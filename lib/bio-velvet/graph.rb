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
              log = Bio::Log::LoggerPlus['bio-velvet']
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
        attr_accessor :node_id, :coverages, :ends_of_kmers_of_node, :ends_of_kmers_of_twin_node

        # For read tracking
        attr_accessor :number_of_short_reads
        # For read tracking - an array of NodedRead objects
        attr_accessor :short_reads

        # The sequence of this node, should a contig be made solely out of this node.
        # The kmer length is that kmer length that was used to create the assembly.
        #
        # If this node has a sequence that is 2 or more less than the hash length, then the
        # sequence of this node requires information outside of this object, and gathering
        # that information is now not implemented here.
        def sequence(kmer_length)
          len = @ends_of_kmers_of_node.length
          if len < kmer_length -1
            raise NotImplementedException, "Attempted to get the sequence of a velvet node that is too short, such that the sequence info is not fully present in the node object"
          end

          # Sequence is the reverse complement of the ends_of_kmers_of_twin_node,
          # Then the ends_of_kmers_of_node after removing the first kmer_length - 1
          # nucleotides
          Bio::Sequence::NA.new(@ends_of_kmers_of_twin_node).reverse_complement.to_s.upcase+
            @ends_of_kmers_of_node[len-kmer_length+1 ... len]
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
      end

      # Tracked read, part of a node
      class NodedRead
        attr_accessor :read_id, :offset_from_start_of_node, :start_coord, :direction
      end
    end
  end
end
