#require 'hopcsv'
require 'bio'
require 'tempfile'

module Bio
  module Velvet
    class NotImplementedException < Exception; end

    # Parser for a velvet assembler's graph file (Graph or LastGraph) output from velvetg
    #
    # The definition of this file is given in the velvet manual, at
    # http://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf
    class Graph
      include Bio::Velvet::Logging

      # Taken directly from the graph, statistics and information about the Graph i.e. from the velvet manual "$NUMBER_OF_NODES $NUMBER_OF_SEQUENCES $HASH_LENGTH"
      attr_accessor :number_of_nodes, :number_of_sequences, :hash_length

      # NodeArray object of all the graph's node objects
      attr_accessor :nodes

      # ArcArray of Arc objects
      attr_accessor :arcs

      def self.log
        self.new.log
      end

      # Parse a graph file from a Graph, Graph2 or LastGraph output file from velvet
      # into a Bio::Velvet::Graph object
      #
      # Options:
      # * :dont_parse_noded_reads: if true, then parsing of the NR section is skipped
      # * :interesting_read_ids: If not nil, is a Set of nodes that we are interested in. Reads
      # not of interest will not be parsed in (the NR part of the velvet LastGraph file). Regardless all
      # nodes and edges are parsed in. Using this options saves both memory and CPU.
      # * :interesting_node_ids: like :interesting_read_ids except it allows targeting of particular nodes
      # rather than particular reads.
      # * :grep_hack: to make the parsing of read associations go even faster, a grep-based, rather
      # hacky method is applied to the graph file, so only NR data of interesting_read_ids is presented
      # to the parser. This can save days of parsing time, but is a bit of a hack and its usage may
      # not be particularly future-proof. The value of this option is the amount of context coming out of grep
      # (the -A and -B flags). Using 500 should probably work for most circumstances, if not an Exception will be raised.
      def self.parse_from_file(path_to_graph_file, options={})
        graph = self.new
        state = :header

        current_node = nil
        graph.nodes = NodeArray.new
        graph.arcs = ArcArray.new
        current_node_direction = nil

        line_number = 0
        Hopcsv.foreach(path_to_graph_file,"\t") do |row|
          line_number += 1

          if state == :header
            raise "parse exception on header line, this line #{line_number}: #{row.inspect}" unless row.length >= 3
            graph.number_of_nodes = row[0].to_i
            graph.number_of_sequences = row[1].to_i
            graph.hash_length = row[2].to_i
            #Not quite sure what the function of the 4th column is
            state = :nodes_0
            log.debug "Now parsing velvet graph nodes" if log.debug?
            next
          end

          if state == :nodes_0
            # NODE $NODE_ID $COV_SHORT1 $O_COV_SHORT1 $COV_SHORT2 $O_COV_SHORT2
            # $ENDS_OF_KMERS_OF_NODE
            # $ENDS_OF_KMERS_OF_TWIN_NODE
            if row[0] == 'NODE'
              raise unless row.length > 2
              current_node = Node.new
              current_node.node_id = row[1].to_i
              current_node.length = row[2].to_i
              current_node.coverages = row[3...row.length].collect{|c| c.to_i}
              current_node.parent_graph = graph
              state = :nodes_1
              raise "Duplicate node name" unless graph.nodes[current_node.node_id].nil?
              graph.nodes[current_node.node_id] = current_node
              next
            else
              state = :arc
              log.debug "Now parsing velvet graph arcs" if log.debug?
              # No next in the loop so that this line gets parsed as an ARC further down the loop
            end
          elsif state == :nodes_1
            # Sometimes nodes can be empty
            row[0] ||= ''
            current_node.ends_of_kmers_of_node = row[0]
            raise "Unexpected nodes_1 type line on line #{line_number}: #{row.inspect}" if row.length != 1
            state = :nodes_2
            next
          elsif state == :nodes_2
            # Sometimes nodes can be empty
            row[0] ||= ''
            raise if row.length != 1
            current_node.ends_of_kmers_of_twin_node = row[0]
            state = :nodes_0
            next
          end

          if state == :arc
            if row[0] == 'ARC'
              # ARC $START_NODE $END_NODE $MULTIPLICITY
              arc = Arc.new
              raise unless row.length == 4
              arc.begin_node_id = row[1].to_i.abs
              arc.end_node_id = row[2].to_i.abs
              arc.multiplicity = row[3].to_i
              arc.begin_node_direction = (row[1].to_i > 0)
              arc.end_node_direction = (row[2].to_i > 0)
              graph.arcs.push arc
              next
            else
              state = :nr
              log.debug "Finished parsing velvet graph arcs. Now parsing the rest of the file" if log.debug?
            end
          end

          if state == :nr
            if row[0] == 'SEQ'
              log.warn "velvet graph parse warning: SEQ lines in the Graph file parsing not implemented yet, tracking of reads now not parsed either"
              break
            end

            # If short reads are tracked, for every node a block of read identifiers:
            # NR $NODE_ID $NUMBER_OF_SHORT_READS
            # $READ_ID $OFFSET_FROM_START_OF_NODE $START_COORD
            # $READ_ID2 etc.
            #p row
            if row[0] == 'NR'
              break if options[:dont_parse_noded_reads] # We are done if NR things aren't parsed
              if options[:grep_hack]
                unless options[:interesting_read_ids] or options[:interesting_node_ids]
                  raise "Programming error using bio-velvet: if :grep_hack is specified, then :interesting_read_ids or :interesting_node_ids must also be"
                end
                apply_grep_hack graph, path_to_graph_file, options[:interesting_read_ids], options[:interesting_node_ids], options[:grep_hack]
                break #no more parsing is required
              else
                raise unless row.length == 3
                node_pm = row[1].to_i
                current_node_direction = node_pm > 0
                current_node = graph.nodes[node_pm.abs]
                current_node.number_of_short_reads ||= 0
                current_node.number_of_short_reads += row[2].to_i
                next
              end
            else
              raise unless row.length == 3
              read_id = row[0].to_i
              if (options[:interesting_node_ids] and !options[:interesting_node_ids].include?(current_node.node_id)) or
                (options[:interesting_read_ids] and !options[:interesting_read_ids].include?(read_id))
                # We have come across an uninteresting read. Ignore it.
                next
              end
              nr = NodedRead.new
              nr.read_id = read_id
              nr.offset_from_start_of_node = row[1].to_i
              nr.start_coord = row[2].to_i
              nr.direction = current_node_direction
              current_node.short_reads ||= []
              current_node.short_reads.push nr
              next
            end
          end
        end
        log.debug "Finished parsing velvet graph file" if log.debug?

        return graph
      end

      # Return an array of Arc objects between two nodes (specified by integer IDs),
      # or an empty array if none exists. There is four possible arcs between
      # two nodes, connecting their beginnings and ends
      def get_arcs_by_node_id(node_id1, node_id2)
        @arcs.get_arcs_by_node_id(node_id1, node_id2)
      end

      # Return an array of Arc objects between two nodes (specified by node objects),
      # or an empty array if none exists. There is four possible arcs between
      # two nodes, connecting their beginnings and ends
      def get_arcs_by_node(node1, node2)
        @arcs.get_arcs_by_node_id(node1.node_id, node2.node_id)
      end

      # Return the adjacent nodes in the graph that connect to the end of a node
      def neighbours_off_end(node)
        # Find all arcs that include this node in the right place
        passable_nodes = []
        @arcs.get_arcs_by_node_id(node.node_id).each do |arc|
          if arc.begin_node_id == node.node_id and arc.begin_node_direction
            # The most intuitive case
            passable_nodes.push nodes[arc.end_node_id]
          elsif arc.end_node_id == node.node_id and !arc.end_node_direction
            passable_nodes.push nodes[arc.begin_node_id]
          end
        end
        return passable_nodes
      end

      # Return the adjacent nodes in the graph that connect to the end of a node
      def neighbours_into_start(node)
        # Find all arcs that include this node in the right place
        passable_nodes = []
        @arcs.get_arcs_by_node_id(node.node_id).each do |arc|
          if arc.end_node_id == node.node_id and arc.end_node_direction
            passable_nodes.push nodes[arc.begin_node_id]
          elsif arc.begin_node_id == node.node_id and !arc.begin_node_direction
            passable_nodes.push nodes[arc.end_node_id]
          end
        end
        return passable_nodes
      end


      # Deletes nodes and associated arcs from the graph if the block passed
      # evaluates to true (as in Array#delete_if). Actually the associated arcs
      # are deleted first, and then the node, so that the graph remains sane at all
      # times - there is never dangling arcs, as such.
      #
      # Returns a [deleted_nodes, deleted_arcs] tuple, which are both enumerables,
      # each in no particular order.
      def delete_nodes_if(&block)
        deleted_nodes = []
        deleted_arcs = []
        nodes.each do |node|
          if yield(node)
            deleted_nodes.push node

            # delete associated arcs
            arcs_to_del = @arcs.get_arcs_by_node_id(node.node_id)
            deleted_arcs.push arcs_to_del
            arcs_to_del.each do |arc|
              @arcs.delete arc
            end

            # delete the arc itself
            nodes.delete node
          end
        end
        return deleted_nodes, deleted_arcs.flatten
      end

      # Add more noded reads to this already parsed graph. There is
      # no gaurantee that old NodedRead information is preserved, or removed.
      #
      # Options:
      # * :interesting_read_ids: If not nil, is a Set of nodes that we are interested in. Reads
      # not of interest will not be parsed in (the NR part of the velvet LastGraph file). Regardless all
      # nodes and edges are parsed in. Using this options saves both memory and CPU.
      # * :interesting_node_ids: like :interesting_read_ids except it allows targeting of particular nodes
      # rather than particular reads.
      # * :grep_hack: to make the parsing of read associations go even faster, a grep-based, rather
      # hacky method is applied to the graph file, so only NR data of interesting_read_ids is presented
      # to the parser. This can save days of parsing time, but is a bit of a hack and its usage may
      # not be particularly future-proof. The value of this option is the amount of context coming out of grep
      # (the -A and -B flags). Using 500 should probably work for most circumstances, if not an Exception will be raised.
      def parse_additional_noded_reads(path_to_graph_file, options)
        grep_context = options[:grep_hack]
        if grep_context.nil?
          raise "Calling Graph#parse_additional_noded_reads without specifying :grep_hack is currently not implemented"
        end
        self.class.apply_grep_hack(self,
          path_to_graph_file,
          options[:interesting_read_ids],
          options[:interesting_node_ids],
          grep_context
          )
      end





      # A container class for a list of Node objects. Can index with 1-offset
      # IDs, so that they line up with the identifiers in velvet Graph files,
      # yet respond sensibly to NodeArray#length, etc.
      class NodeArray
        include Enumerable

        def initialize
          # Internal index is required because when things get deleted stuff changes.
          @internal_structure = {}
        end

        def []=(node_id, value)
          @internal_structure[node_id] = value
        end

        def [](node_id)
          @internal_structure[node_id]
        end

        def delete(node)
          @internal_structure.delete node.node_id
        end

        def length
          @internal_structure.length
        end

        def each(&block)
          @internal_structure.each do |internal_id, node|
            block.yield node
          end
        end
      end

      class ArcArray
        include Enumerable

        def initialize
          # Internal structure is hash of [node_id1, node_id2] => Array of arcs
          @internal_structure = {}
          @node_to_keys = {}
        end

        def push(arc)
          key = [arc.begin_node_id, arc.end_node_id].sort
          @internal_structure[key] ||= []
          @internal_structure[key].push arc
          @node_to_keys[arc.begin_node_id] ||= []
          @node_to_keys[arc.begin_node_id].push key
          unless arc.begin_node_id == arc.end_node_id
            @node_to_keys[arc.end_node_id] ||= []
            @node_to_keys[arc.end_node_id].push key
          end
        end

        # Return all arcs into or out of the given node_id, or
        def get_arcs_by_node_id(node_id1, node_id2=nil)
          if node_id2.nil?
            next_keys = @node_to_keys[node_id1]
            return [] if next_keys.nil?
            next_keys.uniq.collect do |key|
              @internal_structure[key]
            end.flatten
          else
            to_return = @internal_structure[[node_id1, node_id2].sort]
            if to_return.nil?
              return []
            else
              return to_return
            end
          end
        end

        def delete(arc)
          key = [arc.begin_node_id, arc.end_node_id].sort
          @internal_structure[key].delete arc
          # If there is no other arcs with this same key, clean up more
          if @internal_structure[key].empty?
            @internal_structure.delete key
            @node_to_keys[key[0]].delete key
            @node_to_keys[key[1]].delete key
            @node_to_keys[key[0]] = nil if @node_to_keys[key[0]].nil? or @node_to_keys[key[0]].empty?
            @node_to_keys[key[1]] = nil if @node_to_keys[key[1]].nil? or @node_to_keys[key[1]].empty?
          end
        end

        def length
          @internal_structure.values.flatten.length
        end

        def each(&block)
          @internal_structure.each do |internal_id, arcs|
            arcs.each do |arc|
              block.yield arc
            end
          end
        end
      end

      class Node
        include Bio::Velvet::Logging

        attr_accessor :node_id, :coverages, :ends_of_kmers_of_node, :ends_of_kmers_of_twin_node

        # For read tracking
        attr_accessor :number_of_short_reads
        # For read tracking - an array of NodedRead objects
        attr_accessor :short_reads

        # Graph to which this node belongs
        attr_accessor :parent_graph

        # Number of nucleotides in this node if a contig was made from this contig alone
        attr_accessor :length

        # The sequence of this node, should a contig be made solely out of this node.
        # The kmer length is that kmer length that was used to create the assembly.
        #
        # If this node has a sequence that is 2 or more less than the hash length, then the
        # sequence of this node requires information outside of this object, and gathering
        # that information is not implemented here.
        def sequence
          if !sequence?
            raise NotImplementedException, "Attempted to get the sequence of a velvet node that is too short, such that the sequence info is not fully present in the node object"
          end
          kmer_length = @parent_graph.hash_length

           # Sequence is the reverse complement of the ends_of_kmers_of_twin_node,
           # Then the ends_of_kmers_of_node after removing the first kmer_length - 1
           # nucleotides
           length_to_get_from_fwd = corresponding_contig_length - @ends_of_kmers_of_twin_node.length
           fwd_length = @ends_of_kmers_of_node.length
           raise "Programming error" if length_to_get_from_fwd > fwd_length
           revcom(@ends_of_kmers_of_twin_node)+
             @ends_of_kmers_of_node[-length_to_get_from_fwd...fwd_length]
        end

        # Number of nucleotides in this node if this contig length is being added to
        # another node's length (nodes overlap)
        def length_alone
          @ends_of_kmers_of_node.length
        end

        # The common length of [ends_of_kmers_of_node and :ends_of_kmers_of_twin_node]
        # is equal to the length of the corresponding contig minus k − 1.
        #
        # This method returns that corresponding contig's length
        def corresponding_contig_length
          @ends_of_kmers_of_node.length+@parent_graph.hash_length-1
        end

        # Is it possible to extract the sequence of this node? I.e. is it long enough?
        def sequence?
          kmer_length = @parent_graph.hash_length
          if kmer_length -1 > @ends_of_kmers_of_node.length
            return false
          else
            return true
          end
        end

        # The reverse complement of this node's sequence
        def reverse_sequence
          revcom(sequence)
        end

        # Number of nucleotides in this node if this contig length is being added to
        # another node's length (nodes overlap)
        def length_alone
          @ends_of_kmers_of_node.length
        end

        def to_s
          "Node #{@node_id}: #{@ends_of_kmers_of_node} / #{@ends_of_kmers_of_twin_node}"
        end

        def inspect
          to_s
        end

        # Return the sum of all coverage columns, divided by the length of the node,
        # or nil if this node has no coverage
        def coverage
          return nil if length == 0

          coverage = 0
          coverages.each_with_index do |cov, i|
            # Only take the 0th, 2nd, 4th, etc, don't want the O_cov things
            coverage += cov if i.modulo(2) == 0
          end
          return coverage.to_f / length
        end

        private
        def revcom(seq)
          Bio::Sequence::NA.new(seq).reverse_complement.to_s.upcase
        end
      end

      class Arc
        attr_accessor :begin_node_id, :end_node_id, :multiplicity

        # true for forwards direction, false for reverse
        attr_accessor :begin_node_direction, :end_node_direction

        def directions_opposing?
          if (@begin_node_direction == true and @end_node_direction == false) or
            (@begin_node_direction == false and @end_node_direction == true)
            return true
          elsif [true,false].include?(@begin_node_direction) and [true,false].include?(@end_node_direction)
            return false
          else
            raise Exception, "Node directions not set! Cannot tell whether directions are opposing"
          end
        end

        def begin_node_forward?
          @begin_node_direction
        end

        def end_node_forward?
          @end_node_forward
        end

        # Returns true if this arc connects the end of the first node
        # to the start of the second node, else false
        def connects_end_to_beginning?(first_node_id, second_node_id)
          # ARC $START_NODE $END_NODE $MULTIPLICITY
          #Note: this one line implicitly represents an arc from node A to B and
          #another, with same multiplicity, from -B to -A.
          (first_node_id == @begin_node_id and second_node_id == @end_node_id and
            @begin_node_direction == true and @end_node_direction == true) or
            (first_node_id == @end_node_id and second_node_id = @begin_node_id and
            @begin_node_direction == false and @end_node_direction == false)
        end

        # Returns true if this arc connects the end of the first node
        # to the end of the second node, else false
        def connects_end_to_end?(first_node_id, second_node_id)
          (first_node_id == @begin_node_id and second_node_id == @end_node_id and
            @begin_node_direction == true and @end_node_direction == false) or
            (first_node_id == @end_node_id and second_node_id = @begin_node_id and
            @begin_node_direction == true and @end_node_direction == false)
        end

        # Returns true if this arc connects the start of the first node
        # to the start of the second node, else false
        def connects_beginning_to_beginning?(first_node_id, second_node_id)
          (first_node_id == @begin_node_id and second_node_id == @end_node_id and
            @begin_node_direction == false and @end_node_direction == true) or
            (first_node_id == @end_node_id and second_node_id = @begin_node_id and
            @begin_node_direction == false and @end_node_direction == true)
        end

        # Returns true if this arc connects the start of the first node
        # to the start of the second node, else false
        def connects_beginning_to_end?(first_node_id, second_node_id)
          (first_node_id == @begin_node_id and second_node_id == @end_node_id and
            @begin_node_direction == false and @end_node_direction == false) or
            (first_node_id == @end_node_id and second_node_id = @begin_node_id and
            @begin_node_direction == true and @end_node_direction == true)
        end

        # Return true if this arc connects the beginning of the node,
        # else false
        def connects_to_beginning?(node_id)
          (node_id == @begin_node_id and !@begin_node_direction) or
          (node_id == @end_node_id and @end_node_direction)
        end

        # Return true if this arc connects the end of the node,
        # else false
        def connects_to_end?(node_id)
          (node_id == @begin_node_id and @begin_node_direction) or
          (node_id == @end_node_id and !@end_node_direction)
        end

        def to_s
          str = ''
          str += '-' if @begin_node_direction == false
          str += @begin_node_id.to_s
          str += ' '
          str += '-' if @end_node_direction == false
          str += @end_node_id.to_s
          str += ' '
          str += @multiplicity.to_s
          str
        end
      end

      # Tracked read, part of a node
      class NodedRead
        attr_accessor :read_id, :offset_from_start_of_node, :start_coord, :direction
      end

      private
      def self.apply_grep_hack(graph, path_to_graph_file, interesting_read_ids, interesting_node_ids, grep_context)
        interesting_read_ids ||= []
        interesting_node_ids ||= []
        if interesting_read_ids.empty? and interesting_node_ids.empty?
          log.debug "Nothing to grep for in grep hack" if log.debug?
          return
        end

        Tempfile.open('grep_v_hack') do |tempfile|
          # Create a file to pass to grep -f
          unless interesting_read_ids.nil?
            interesting_read_ids.each do |read_id|
              tempfile.puts "^#{read_id}\t"
            end
          end
          unless interesting_node_ids.nil?
            interesting_node_ids.each do |node_id|
              tempfile.puts "^NR\t#{node_id}\t"
              tempfile.puts "^NR\t-#{node_id}\t"
            end
          end
          tempfile.close

          cmd = "grep -B #{grep_context.inspect} -A #{grep_context.inspect} -f #{tempfile.path} #{path_to_graph_file.inspect}"
          # TODO: make this call more robust
          # grep_result = Bio::Commandeer.run cmd
          s, grep_result, stderr = systemu cmd

          # Parse the grepped out results
          current_node = nil
          current_node_direction = nil
          in_nr_section = false
          grep_result.each_line do |line|
            row = line.split("\t")
            if in_nr_section == false
              # If there is a lot of context then the context includes ARC definitions etc. Skip past this.
              if row[0] == 'NR'
                in_nr_section = true
              elsif row[0] == '--'
                raise "Parsing exception - grep hack too hacky. Sorry. Try modifying the code to increase the default amount of context grep is giving"
              else
                next #skip to next line, waiting to ge into NR section
              end
            end

            if line == "--\n" #the break introduced by grep
              # If we encounter a grep break, but haven't assigned any nodes, then that's not good enough
              if current_node.nil?
                raise "Parsing exception - grep hack too hacky. Sorry. Try modifying the code to increase the default amount of context grep is giving"
              end
              # reset the parsing situation
              current_node = nil
            elsif row[0] == 'NR'
              raise unless row.length == 3
              node_pm = row[1].to_i
              current_node_direction = node_pm > 0
              current_node = graph.nodes[node_pm.abs]
              current_node.number_of_short_reads ||= 0
              current_node.number_of_short_reads += row[2].to_i
              next
            else
              raise unless row.length == 3
              read_id = row[0].to_i
              if (current_node.nil? or !interesting_node_ids.include?(current_node.node_id)) and
                !interesting_read_ids.include?(read_id)
                # We have come across an uninteresting read. Ignore it.
                next
              end
              if current_node.nil?
                # Came across a high coverage node, and grep isn't giving enough context. Hopefully this won't happen much
                # particularly if the reads you are interested in are given to velvet first
                raise "Parsing exception - grep hack too hacky. Sorry. Try modifying the code to increase the default amount of context grep is giving"
              end
              nr = NodedRead.new
              nr.read_id = read_id
              nr.offset_from_start_of_node = row[1].to_i
              nr.start_coord = row[2].to_i
              nr.direction = current_node_direction
              current_node.short_reads ||= []
              current_node.short_reads.push nr
              next
            end
          end

          if current_node.nil?
            raise "Parsing exception - grep hack too hacky. Sorry. Try modifying the code to increase the default amount of context grep is giving"
          end
        end
      end
    end
  end
end
