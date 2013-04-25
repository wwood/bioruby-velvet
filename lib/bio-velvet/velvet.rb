require 'csv'

module Bio
  module Velvet
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

      # Parse a graph file from a Graph or LastGraph output file from velvet
      def self.parse_from_file(path_to_graph_file)
        graph = self.new
        state = :header

        current_node = nil
        graph.nodes = []
        graph.arcs = []

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
              # No next in the loop so that this line gets parsed later on
            end
          elsif state == :nodes_1
            raise if row.length != 1
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
            else
              # There is more possible, but that's enough for my purposes for the moment
              break
            end
          end
        end

        return graph
      end

      class Node
        attr_accessor :node_id, :coverages, :ends_of_kmers_of_node, :ends_of_kmers_of_twin_node
      end

      class Arc
        attr_accessor :begin_node_id, :end_node_id, :multiplicity

        # true for forwards direction, false for reverse
        attr_accessor :begin_node_direction, :end_node_direction
      end
    end
  end
end
