require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe "BioVelvet" do
  it "should be able to parse a graph" do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_reads_assembly','Graph')
    graph.should be_kind_of(Bio::Velvet::Graph)

    graph.number_of_nodes.should eq(967)
    graph.number_of_sequences.should eq(50000)
    graph.hash_length.should eq(31)

    graph.nodes[1].should be_kind_of(Bio::Velvet::Graph::Node)
    graph.nodes.length.should eq(967+1) #+1 because it is stored in an array object
    graph.nodes[1].node_id.should eq(1)
    graph.nodes[3].coverages.should eq([3,236,205,0,0])
    graph.nodes[3].ends_of_kmers_of_node.should eq('TTG')
    graph.nodes[3].ends_of_kmers_of_twin_node.should eq('ACA')

    graph.arcs.length.should eq(563)
    graph.arcs[0].begin_node_id.should eq(2)
    graph.arcs[0].end_node_id.should eq(712)
    graph.arcs[0].begin_node_direction.should eq(true)
    graph.arcs[0].end_node_direction.should eq(false)
    graph.arcs[0].multiplicity.should eq(1)

    graph.arcs[2].begin_node_direction.should eq(false)
  end
end
