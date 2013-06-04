require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'

describe "BioVelvet" do
  it "should be able to parse a graph" do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_reads_assembly','Graph')
    graph.should be_kind_of(Bio::Velvet::Graph)

    graph.number_of_nodes.should eq(967)
    graph.number_of_sequences.should eq(50000)
    graph.hash_length.should eq(31)

    graph.nodes[1].should be_kind_of(Bio::Velvet::Graph::Node)
    graph.nodes.length.should eq(967)
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

  it "should be able to parse a graph that has read tracking" do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_reads_assembly_read_tracking','Graph2')
    graph.should be_kind_of(Bio::Velvet::Graph)

    graph.number_of_nodes.should eq(967)
    graph.number_of_sequences.should eq(50000)
    graph.hash_length.should eq(31)

    graph.nodes[1].should be_kind_of(Bio::Velvet::Graph::Node)
    graph.nodes.length.should eq(967)
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

    # NR	-967	1
    # 49982	0	0
    # === later
    # NR	967	1
    # 49981	0	0
    node = graph.nodes[967]
    node.short_reads.nil?.should eq(false)
    node.short_reads.length.should eq(2), node.inspect
    node.short_reads[0].read_id.should eq(49982)
    node.short_reads[0].offset_from_start_of_node.should eq(0)
    node.short_reads[0].start_coord.should eq(0)
    node.short_reads[0].direction.should eq(false)
    node.short_reads[1].read_id.should eq(49981)
    node.short_reads[1].offset_from_start_of_node.should eq(0)
    node.short_reads[1].start_coord.should eq(0)
    node.short_reads[1].direction.should eq(true)

    # NR	-951	2
    #47210	0	0
    #47223	41	0
    # ====later
    # NR	951	2
    # 47209	54	0
    # 47224	0	0
    node = graph.nodes[951]
    node.short_reads.length.should eq(4)
    node.number_of_short_reads.should eq(4)
    node.short_reads[1].offset_from_start_of_node.should eq(41)

    # grep -A 50000 ^NR velvet_test_reads_assembly_read_tracking/Graph2 |grep -v NR |wc -l
    # 40327
    graph.nodes.collect{|n| n.short_reads.nil? ? 0 : n.short_reads.length}.reduce(:+).should eq(40327)
  end

  it 'should have a functioning NodeArray class' do
    na = Bio::Velvet::Graph::NodeArray.new
    na.length.should eq(0)
    node = Bio::Velvet::Graph::Node.new
    na[1] = node
    na.length.should eq(1)
    na[1].should eq(node)
  end

  it "arcs should directions_opposing?" do
    arc = Bio::Velvet::Graph::Arc.new
    arc.begin_node_direction = true
    arc.end_node_direction = true
    arc.directions_opposing?.should eq(false)

    arc.begin_node_direction = true
    arc.end_node_direction = false
    arc.directions_opposing?.should eq(true)

    arc.begin_node_direction = false
    arc.end_node_direction = false
    arc.directions_opposing?.should eq(false)
  end

  it "nodes should correctly respond to #sequence" do
    node = Bio::Velvet::Graph::Node.new
    node.ends_of_kmers_of_node = 'AATCAAACTATAAAGTTTTAGAAGATAAAGTAACAACTTATACATGGGGATTCGGAGTTAAAAAAGTAGATTCAGAAAATATTTCAATAGATCTTGCAGGCGCAGCTTTTTCTGTTAGGGATAAAAATGGTAATGTAATTGGTAAATATACGTATGATTCTACTGGAAATGTGGTTTTATTAAAAGGAAAGGGTGTAACTGATAAAAATGGACGAGTTATATTTACTGGTTTAAAAGAAGGAGATTACTTTATAAAAGAAGAAAAAGCTCCTAAAGGGTATAGCCTTTTAAAAGAACCAGTAAAAGTTACTATAACAGCTCAAAAAGATGATAATGGAGAGTATACTGGTCAAGCAACTATATCTGTAACTAATGGCAATGAAGCTGGAAGTATAATAAATAATATTACTATGAATGATGGCAATGTATTATTTAATGTACAAATTAAAAACTATGCTGGTATTTCACTTCCAGGTACAGG'
    node.ends_of_kmers_of_twin_node = 'TTTTTAATTTGTACATTAAATAATACATTGCCATCATTCATAGTAATATTATTTATTATACTTCCAGCTTCATTGCCATTAGTTACAGATATAGTTGCTTGACCAGTATACTCTCCATTATCATCTTTTTGAGCTGTTATAGTAACTTTTACTGGTTCTTTTAAAAGGCTATACCCTTTAGGAGCTTTTTCTTCTTTTATAAAGTAATCTCCTTCTTTTAAACCAGTAAATATAACTCGTCCATTTTTATCAGTTACACCCTTTCCTTTTAATAAAACCACATTTCCAGTAGAATCATACGTATATTTACCAATTACATTACCATTTTTATCCCTAACAGAAAAAGCTGCGCCTGCAAGATCTATTGAAATATTTTCTGAATCTACTTTTTTAACTCCGAATCCCCATGTATAAGTTGTTACTTTATCTTCTAAAACTTTATAGTTTGATTCTAAATCGTGATCTTTGGTAGAGATAAGTG'
    node.sequence(31).should eq('CACTTATCTCTACCAAAGATCACGATTTAGAATCAAACTATAAAGTTTTAGAAGATAAAG
TAACAACTTATACATGGGGATTCGGAGTTAAAAAAGTAGATTCAGAAAATATTTCAATAG
ATCTTGCAGGCGCAGCTTTTTCTGTTAGGGATAAAAATGGTAATGTAATTGGTAAATATA
CGTATGATTCTACTGGAAATGTGGTTTTATTAAAAGGAAAGGGTGTAACTGATAAAAATG
GACGAGTTATATTTACTGGTTTAAAAGAAGGAGATTACTTTATAAAAGAAGAAAAAGCTC
CTAAAGGGTATAGCCTTTTAAAAGAACCAGTAAAAGTTACTATAACAGCTCAAAAAGATG
ATAATGGAGAGTATACTGGTCAAGCAACTATATCTGTAACTAATGGCAATGAAGCTGGAA
GTATAATAAATAATATTACTATGAATGATGGCAATGTATTATTTAATGTACAAATTAAAA
ACTATGCTGGTATTTCACTTCCAGGTACAGG'.gsub(/\n/,''))

    #having problems am I? It seemed so.
    node.ends_of_kmers_of_node = 'GTTTAAAAGAAGGAGATTACTTTATAAAA'
    node.ends_of_kmers_of_twin_node = 'AGTAAATATAACTCGTCCATTTTTATCAG'
    #  lambda {walker.trail_sequence(graph, [nodes[2],nodes[4],nodes[3]]).should raise_error(Bio::AssemblyGraphAlgorithms::GraphWalkingException)}
    lambda {node.sequence(31).should raise_error(Bio::Velvet::NotImplementedException)}
  end
end
