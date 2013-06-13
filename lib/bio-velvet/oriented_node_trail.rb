module Bio
  module Velvet
    class Graph
      # An ordered list of nodes, each with an orientation along that trail
      class OrientedNodeTrail
        include Enumerable

        START_IS_FIRST = :start_is_first
        END_IS_FIRST = :end_is_first

        def initialize
          @trail = []
        end

        # Add a node to the trail. start_or_end is either
        # OrientedNodeTrail::START_IS_FIRST or OrientedNodeTrail::END_IS_FIRST
        def add_node(node, start_or_end)
          oriented = OrientedNode.new
          oriented.node = node
          oriented.first_side = start_or_end
          trail.push oriented
        end

        def each(&block)
          @trail.each(&block)
        end

        # Return the sequence of the entire trail, or an empty string if there is no
        # nodes in the trail.
        def sequence
          return '' if @trail.empty?
          initial = nil
          if @trail[0].starts_at_start?
            initial = node.sequence
          else
            initial = node.reverse_sequence
          end

          return @trail.reduce(@trail[0].ends_of_kmers_of_twin_node) do |prev, oriented|
            if @trail.starts_at_start?
              prev + oriented.ends_of_kmers_of_node
            else
              prev + oriented.ends_of_kmers_of_twin_node
            end
          end
        end

        class OrientedNode
          attr_accessor :node, :first_side

          def starts_at_start?
            @first_side == OrientedNodeTrail::START_IS_FIRST
          end

          def starts_at_end?
            @first_side == OrientedNodeTrail::END_IS_FIRST
          end
        end
      end
    end
  end
end
