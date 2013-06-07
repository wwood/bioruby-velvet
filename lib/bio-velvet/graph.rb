require 'csv'
require 'bio'

module Bio
  module Velvet
    class NotImplementedException < Exception; end

    # Parser for a velvet assembler's graph file (Graph or LastGraph) output from velvetg
    #
    # The definition of this file is given in the velvet manual, at
    # http://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf
    class Graph
      include Bio::Velvet::Logging

      # $NUMBER_OF_NODES $NUMBER_OF_SEQUENCES $HASH_LENGTH
      attr_accessor :number_of_nodes, :number_of_sequences, :hash_length

      # Hash of node identifying integers to Node objects
      attr_accessor :nodes

      # Array of Arc objects
      attr_accessor :arcs

      # Parse a graph file from a Graph, Graph2 or LastGraph output file from velvet
      # into a Bio::Velvet::Graph object
      def self.parse_from_file(path_to_graph_file)
        graph = self.new
        state = :header

        current_node = nil
        graph.nodes = NodeArray.new
        graph.arcs = []
        current_node_direction = nil

        CSV.foreach(path_to_graph_file, :col_sep => "\t") do |row|
          if state == :header
            raise "parse exception on header line, this line: #{row.inspect}" unless row.length >= 3
            graph.number_of_nodes = row[0].to_i
            graph.number_of_sequences = row[1].to_i
            graph.hash_length = row[2].to_i
            state = :nodes_0
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
              current_node.coverages = row[2...row.length].collect{|c| c.to_i}
              current_node.parent_graph = graph
              state = :nodes_1
              raise "Duplicate node name" unless graph.nodes[current_node.node_id].nil?
              graph.nodes[current_node.node_id] = current_node
              next
            else
              state = :arc
              # No next in the loop so that this line gets parsed as an ARC further down the loop
            end
          elsif state == :nodes_1
            raise "Unexpected nodes_1 type line: #{row.inspect}" if row.length != 1
            current_node.ends_of_kmers_of_node = row[0]
            state = :nodes_2
            next
          elsif state == :nodes_2
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
              raise unless row.length == 3
              node_pm = row[1].to_i
              current_node_direction = node_pm > 0
              current_node = graph.nodes[node_pm.abs]
              current_node.number_of_short_reads ||= 0
              current_node.number_of_short_reads += row[2].to_i
              next
            else
              raise unless row.length == 3
              nr = NodedRead.new
              nr.read_id = row[0].to_i
              nr.offset_from_start_of_node = row[1].to_i
              nr.start_coord = row[2].to_i
              nr.direction = current_node_direction
              current_node.short_reads ||= []
              current_node.short_reads.push nr
              next
            end
          end
        end

        return graph
      end

      # Return an array of Arc objects between two nodes (specified by integer IDs),
      # or an empty array if none exists. There is four possible arcs between
      # two nodes, connecting their beginnings and ends
      def get_arcs_by_node_id(node1, node2)
        arcs = []
        @arcs.each do |arc|
          if (arc.begin_node_id == node1 and arc.end_node_id == node2) or
            (arc.begin_node_id == node2 and arc.end_node_id == node1)
            arcs.push arc
          end
        end
        return arcs
      end

      # Return an array of Arc objects between two nodes (specified by node objects),
      # or an empty array if none exists. There is four possible arcs between
      # two nodes, connecting their beginnings and ends
      def get_arcs_by_node(node1, node2)
        arcs = []
        @arcs.each do |arc|
          if (arc.begin_node_id == node1.node_id and arc.end_node_id == node2.node_id) or
            (arc.begin_node_id == node2.node_id and arc.end_node_id == node1.node_id)
            arcs.push arc
          end
        end
        return arcs
      end

      # Return the adjacent nodes in the graph that connect to the end of a node
      def neighbours_off_end(node)
        # Find all arcs that include this node in the right place
        passable_arcs = []
        arcs.each do |arc|
          if arc.begin_node_id == node.node_id and arc.begin_node_direction
            # The most intuitive case
            passable_arcs.push nodes[arc.end_node_id]
          elsif arc.end_node_id == node.node_id and !arc.begin_node_direction
            passable_arcs.push nodes[arc.begin_node_id]
          end
        end
        return passable_arcs
      end

      # Return the adjacent nodes in the graph that connect to the end of a node
      def neighbours_into_start(node)
        # Find all arcs that include this node in the right place
        passable_arcs = []
        arcs.each do |arc|
          if arc.end_node_id == node.node_id and arc.begin_node_direction
            passable_arcs.push nodes[arc.begin_node_id]
          elsif arc.begin_node_id == node.node_id and !arc.begin_node_direction
            passable_arcs.push nodes[arc.end_node_id]
          end
        end
        return passable_arcs
      end





      # A container class for a list of Node objects. Can index with 1-offset
      # IDs, so that they line up with the identifiers in velvet Graph files,
      # yet respond sensibly to NodeArray#length, etc.
      class NodeArray < Array
        def []=(node_id, value)
          raise if node_id < 1
          super(node_id-1, value)
        end

        def [](index)
          super(index-1)
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

        # The sequence of this node, should a contig be made solely out of this node.
        # The kmer length is that kmer length that was used to create the assembly.
        #
        # If this node has a sequence that is 2 or more less than the hash length, then the
        # sequence of this node requires information outside of this object, and gathering
        # that information is now not implemented here.
        def sequence
          kmer_length = @parent_graph.hash_length
          len = @ends_of_kmers_of_node.length
          if len < kmer_length -1 #if the node sequence information cannot be derived from the data in this current object
            # Find an adjacent node
            sequence_length_to_get = 2*kmer_length-len-3
            log.debug "this node has to go looking for missing sequence, need #{sequence_length_to_get} added to either side" if log
            neighbour = @parent_graph.neighbours_into_start(self).max{|a,b| a.ends_of_kmers_of_node.length}

            if !neighbour.nil? and neighbour.ends_of_kmers_of_node.length >= sequence_length_to_get
              log.debug "Found a neighbour into start of this node: #{neighbour.node_id}" if log and log.debug?
              # There's a node coming into this node's start. The sequence I want is the last (hash_length-1)
              # nucleotides of that sequence.
              arcs = @parent_graph.get_arcs_by_node_id(neighbour.node_id, @node_id)
              # There must be at least 1 arc here, otherwise neighbours_into_start won't have returned anything
              arc = arcs[0]
              if arc.connects_end_to_beginning?(neighbour.node_id, @node_id)
                neighbour_seq = neighbour.ends_of_kmers_of_node
                return neighbour_seq[neighbour_seq.length-sequence_length_to_get...neighbour_seq.length]+@ends_of_kmers_of_node
              elsif arc.connects_beginning_to_beginning?(neighbour.node_id, @node_id)
                neighbour_seq = revcom(neighbour.ends_of_kmers_of_twin_node[0...sequence_length_to_get])
                return neighbour_seq+@ends_of_kmers_of_node
              else
                raise "Programming error, or unexpected/malformed velvet graph format file. Node id #{@node_id}"
              end
            else
              neighbour = @parent_graph.neighbours_off_end(self).max{|a,b| a.ends_of_kmers_of_node.length}
              if !neighbour.nil? and neighbour.ends_of_kmers_of_node.length >= sequence_length_to_get
                log.debug "Found a neighbour off end of this node: #{neighbour.node_id}" if log and log.debug?
                arcs = @parent_graph.get_arcs_by_node_id(neighbour.node_id, @node_id)
                # There must be at least 1 arc here, otherwise neighbours_off_end won't have returned anything
                arc = arcs[0]
                if arc.connects_end_to_end?(neighbour.node_id, @node_id)
                  log.debug "Adding node end to end: #{neighbour}"
                  # Add the last bit of the fwd seq of the neighbour (revcom'd)
                  # To the end of the current twin node, and then revcom the
                  # entire thing
                  neighbour_seq = neighbour.ends_of_kmers_of_node
                  neighbour_seq_revcom = revcom(neighbour_seq)[0...sequence_length_to_get]
                  return revcom(@ends_of_kmers_of_twin_node)+neighbour_seq_revcom
                elsif arc.connects_beginning_to_end?(neighbour.node_id, @node_id)
                  # Add the first bit of the twin node to the end of the current twin node
                  # and then revcom the entire thing
                  log.debug "Adding node beginning to end: #{neighbour}"
                  neighbour_seq = neighbour.ends_of_kmers_of_twin_node
                  return revcom(@ends_of_kmers_of_twin_node+neighbour_seq[0...sequence_length_to_get])
                else
                  raise "Programming error, or unexpected/malformed velvet graph format file. Node id #{@node_id}"
                end
              else
                raise NotImplementedException, "Attempting to get the sequence of a short node whose neighbours are also short. This could be implemented in the code, but I'm being lazy for now, hoping it doesn't happen. Node id #{@node_id}"
              end
            end
          else
            # There is sufficient local information for the sequence to be found.
            #
            # Sequence is the reverse complement of the ends_of_kmers_of_twin_node,
            # Then the ends_of_kmers_of_node after removing the first kmer_length - 1
            # nucleotides
            return Bio::Sequence::NA.new(@ends_of_kmers_of_twin_node).reverse_complement.to_s.upcase+
              @ends_of_kmers_of_node[len-kmer_length+1 ... len]
          end
        end

        # Number of nucleotides in this node if a contig was made from this contig alone
        def length
          @parent_graph.hash_length-1+@ends_of_kmers_of_node.length
        end

        # Number of nucleotides in this node if this contig length is being added to
        # another node's length (nodes overlap)
        def length_alone
          @ends_of_kmers_of_node.length
        end

        def to_s
          "Node #{@node_id}: #{@ends_of_kmers_of_node} / #{@ends_of_kmers_of_twin_node}"
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
    end
  end
end
