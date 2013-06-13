require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'

class String
  def revcom
    Bio::Sequence::NA.new(self).reverse_complement.to_s.upcase
  end
end

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
    graph.nodes[3].length.should == 3
    graph.nodes[3].coverages.should eq([236,205])
    graph.nodes[3].ends_of_kmers_of_node.should eq('TTG')
    graph.nodes[3].ends_of_kmers_of_twin_node.should eq('ACA')

    graph.arcs.length.should eq(563)
    arcs = graph.get_arcs_by_node_id(2,712)
    arcs.length.should == 1
    arcs[0].begin_node_id.should eq(2)
    arcs[0].end_node_id.should eq(712)
    arcs[0].begin_node_direction.should eq(true)
    arcs[0].end_node_direction.should eq(false)
    arcs[0].multiplicity.should eq(1)
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
    graph.nodes[3].length.should == 3
    graph.nodes[3].coverages.should eq([236,205])
    graph.nodes[3].ends_of_kmers_of_node.should eq('TTG')
    graph.nodes[3].ends_of_kmers_of_twin_node.should eq('ACA')


    arcs = graph.get_arcs_by_node_id(2,712)
    arcs.length.should == 1
    arcs[0].begin_node_id.should eq(2)
    arcs[0].end_node_id.should eq(712)
    arcs[0].begin_node_direction.should eq(true)
    arcs[0].end_node_direction.should eq(false)
    arcs[0].multiplicity.should eq(1)

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

  it 'should return sets of arcs by id' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_reads_assembly','LastGraph')
    #    ARC     2       -578    1
    #    ARC     2       -473    30
    #    ARC     -2      650     3
    #    ARC     -2      959     24
    #    ARC     3       4       81
    #    ARC     -3      -786    61
    #    ARC     -3      -740    1
    #    ARC     -3      -568    1
    #    ARC     -3      754     6
    # ....
    #Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('info'); log = Bio::Log::LoggerPlus.new('bio-velvet'); Bio::Log::CLI.configure('bio-velvet')

    arcs = graph.get_arcs_by_node_id(2,578)
    arcs.length.should eq(1)
    arcs[0].begin_node_forward?.should eq(true)
    arcs[0].begin_node_id.should == 2
    arcs[0].end_node_id.should == 578

    arcs = graph.get_arcs_by_node_id(578,2)
    arcs.length.should eq(1)
    arcs[0].begin_node_forward?.should eq(true)
    arcs[0].begin_node_id.should == 2
    arcs[0].end_node_id.should == 578

    arcs = graph.get_arcs_by_node_id(2,178)
    arcs.length.should == 0
  end

  it 'should return a set of arcs by node objects' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'velvet_test_reads_assembly','LastGraph')
    #    ARC     2       -578    1
    #    ARC     2       -473    30
    #    ARC     -2      650     3
    #    ARC     -2      959     24
    #    ARC     3       4       81
    #    ARC     -3      -786    61
    #    ARC     -3      -740    1
    #    ARC     -3      -568    1
    #    ARC     -3      754     6
    # ....
    #Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('info'); log = Bio::Log::LoggerPlus.new('bio-velvet'); Bio::Log::CLI.configure('bio-velvet')

    node2 = graph.nodes.select{|n| n.node_id == 2}[0]
    node650 = graph.nodes.select{|n| n.node_id == 650}[0]
    node754 = graph.nodes.select{|n| n.node_id == 754}[0]

    # forward
    arcs = graph.get_arcs_by_node(node2, node650)
    arcs.length.should eq(1)
    arcs[0].begin_node_forward?.should eq(false)
    arcs[0].begin_node_id.should == 2
    arcs[0].end_node_id.should == 650

    #reverse
    arcs = graph.get_arcs_by_node(node650, node2)
    arcs.length.should eq(1)
    arcs[0].begin_node_forward?.should eq(false)
    arcs[0].begin_node_id.should == 2
    arcs[0].end_node_id.should == 650

    # no connection
    arcs = graph.get_arcs_by_node(node2, node754)
    arcs.length.should == 0
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
    graph = Bio::Velvet::Graph.new
    graph.hash_length = 31
    node = Bio::Velvet::Graph::Node.new
    node.parent_graph = graph
    node.ends_of_kmers_of_node = 'AATCAAACTATAAAGTTTTAGAAGATAAAGTAACAACTTATACATGGGGATTCGGAGTTAAAAAAGTAGATTCAGAAAATATTTCAATAGATCTTGCAGGCGCAGCTTTTTCTGTTAGGGATAAAAATGGTAATGTAATTGGTAAATATACGTATGATTCTACTGGAAATGTGGTTTTATTAAAAGGAAAGGGTGTAACTGATAAAAATGGACGAGTTATATTTACTGGTTTAAAAGAAGGAGATTACTTTATAAAAGAAGAAAAAGCTCCTAAAGGGTATAGCCTTTTAAAAGAACCAGTAAAAGTTACTATAACAGCTCAAAAAGATGATAATGGAGAGTATACTGGTCAAGCAACTATATCTGTAACTAATGGCAATGAAGCTGGAAGTATAATAAATAATATTACTATGAATGATGGCAATGTATTATTTAATGTACAAATTAAAAACTATGCTGGTATTTCACTTCCAGGTACAGG'
    node.ends_of_kmers_of_twin_node = 'TTTTTAATTTGTACATTAAATAATACATTGCCATCATTCATAGTAATATTATTTATTATACTTCCAGCTTCATTGCCATTAGTTACAGATATAGTTGCTTGACCAGTATACTCTCCATTATCATCTTTTTGAGCTGTTATAGTAACTTTTACTGGTTCTTTTAAAAGGCTATACCCTTTAGGAGCTTTTTCTTCTTTTATAAAGTAATCTCCTTCTTTTAAACCAGTAAATATAACTCGTCCATTTTTATCAGTTACACCCTTTCCTTTTAATAAAACCACATTTCCAGTAGAATCATACGTATATTTACCAATTACATTACCATTTTTATCCCTAACAGAAAAAGCTGCGCCTGCAAGATCTATTGAAATATTTTCTGAATCTACTTTTTTAACTCCGAATCCCCATGTATAAGTTGTTACTTTATCTTCTAAAACTTTATAGTTTGATTCTAAATCGTGATCTTTGGTAGAGATAAGTG'
    node.sequence.should eq('CACTTATCTCTACCAAAGATCACGATTTAGAATCAAACTATAAAGTTTTAGAAGATAAAG
TAACAACTTATACATGGGGATTCGGAGTTAAAAAAGTAGATTCAGAAAATATTTCAATAG
ATCTTGCAGGCGCAGCTTTTTCTGTTAGGGATAAAAATGGTAATGTAATTGGTAAATATA
CGTATGATTCTACTGGAAATGTGGTTTTATTAAAAGGAAAGGGTGTAACTGATAAAAATG
GACGAGTTATATTTACTGGTTTAAAAGAAGGAGATTACTTTATAAAAGAAGAAAAAGCTC
CTAAAGGGTATAGCCTTTTAAAAGAACCAGTAAAAGTTACTATAACAGCTCAAAAAGATG
ATAATGGAGAGTATACTGGTCAAGCAACTATATCTGTAACTAATGGCAATGAAGCTGGAA
GTATAATAAATAATATTACTATGAATGATGGCAATGTATTATTTAATGTACAAATTAAAA
ACTATGCTGGTATTTCACTTCCAGGTACAGG'.gsub(/\n/,''))
  end

  it 'short nodes should respond to sequence properly' do
    #Bio::Log::CLI.logger('stderr'); Bio::Log::CLI.trace('debug'); log = Bio::Log::LoggerPlus.new('bio-velvet'); Bio::Log::CLI.configure('bio-velvet')
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR,'short_node_LastGraph')
    graph.nodes.length.should == 4
    graph.nodes[2].sequence.should == 'CTGATAAAAATGGACGAGTTATATTTACTG'+'GTTTAAAAGAAGGAGATTACTTTATAAAA'

    # Now test when going from the start cannot be done
    to_delete = nil
    graph.get_arcs_by_node_id(1,2).each{|arc| graph.arcs.delete arc}
    raise unless graph.arcs.length == 3
    graph.nodes[2].sequence.should == 'CTGATAAAAATGGACGAGTTATATTTACTG'+'GTTTAAAAGAAGGAGATTACTTTATAAAA'

    # TODO: only 2 of the 4 if statemenet
  end

  it 'should really correctly respond to #sequence on real data' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'node_sequence','LastGraph')
    graph.should be_kind_of(Bio::Velvet::Graph)

    contig_seq_hash = {}
    Bio::FlatFile.foreach(File.join(TEST_DATA_DIR, 'node_sequence','contigs.fa')) do |seq|
      node_id = seq.definition.match(/^NODE_(\d+)_/)[1]
      contig_seq_hash[node_id] = seq.seq
    end
    graph.nodes.each do |node|
      seq = node.sequence
      exp = contig_seq_hash[node.node_id.to_s]
      seq.should eq(exp), "Node #{node.node_id}. Found\n#{seq}\nExpected\n#{exp}"
    end
  end

  it 'should delete nodes correctly path nodes' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'short_node_LastGraph')
    graph.should be_kind_of(Bio::Velvet::Graph)

    deleted_nodes, deleted_arcs = graph.delete_nodes_if do |node|
      node.node_id == 3
    end
    deleted_nodes.collect{|n| n.node_id}.sort.should == [3]
    graph.nodes.collect{|n| n.node_id}.sort.should == [1,2,4]
    deleted_arcs.collect{|a| [a.begin_node_id, a.end_node_id]}.sort.should == [[1,3],[3,4]]
    graph.arcs.collect{|a| [a.begin_node_id, a.end_node_id]}.sort.should == [[1,2],[2,4]]

    graph.nodes[1].node_id.should == 1
    graph.nodes[3].should == nil
    graph.nodes[4].node_id.should == 4
  end

  it 'should delete nodes correctly cap nodes' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'short_node_LastGraph')
    graph.should be_kind_of(Bio::Velvet::Graph)

    deleted_nodes, deleted_arcs = graph.delete_nodes_if do |node|
      node.node_id == 4
    end
    deleted_nodes.collect{|n| n.node_id}.sort.should == [4]
    graph.nodes.collect{|n| n.node_id}.sort.should == [1,2,3]
    deleted_arcs.collect{|a| [a.begin_node_id, a.end_node_id]}.sort.should == [[2,4],[3,4]]
    graph.arcs.collect{|a| [a.begin_node_id, a.end_node_id]}.sort.should == [[1,2],[1,3]]

    graph.nodes[4].should == nil
    graph.nodes[1].node_id.should == 1
  end

  it 'node#coverage' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'short_node_LastGraph')
    graph.should be_kind_of(Bio::Velvet::Graph)

    graph.nodes[1].coverage.should == 1140.0/228
    graph.nodes[2].coverage.should == 58.0/29
    graph.nodes[3].coverage.should == 114.0/38
    graph.nodes[4].coverage.should == 1120.0/224
  end

  it 'should get neighbours_off_end' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'short_node_LastGraph')
    graph.should be_kind_of(Bio::Velvet::Graph)

    graph.neighbours_off_end(graph.nodes[1]).should == [graph.nodes[2], graph.nodes[3]]
    graph.neighbours_off_end(graph.nodes[2]).should == [graph.nodes[4]]
    graph.neighbours_off_end(graph.nodes[3]).should == [graph.nodes[4]]
    graph.neighbours_off_end(graph.nodes[4]).should == [graph.nodes[2], graph.nodes[3]] #node 4 is revcom
  end

  it 'should get neighbours_into_start' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'short_node_LastGraph')
    graph.should be_kind_of(Bio::Velvet::Graph)

    graph.neighbours_into_start(graph.nodes[1]) == []
    graph.neighbours_into_start(graph.nodes[2]) == [graph.nodes[1]]
    graph.neighbours_into_start(graph.nodes[3]) == [graph.nodes[1]]
    graph.neighbours_into_start(graph.nodes[4]) == []
  end

  it 'should get node sequences properly, when sequences are short and have short links' do
    graph = Bio::Velvet::Graph.parse_from_file File.join(TEST_DATA_DIR, 'short_node_sequence_test_graph')
    graph.should be_kind_of(Bio::Velvet::Graph)

    exp = 'TCCTCAGCATGTTTGTTATTTATGACGTATCAAACAAAACTGAGGACATGTGAAATAACCCCGAATGAGAATATATGCTTTCCTATTGGAACGATTCTTGCTGTAAAAAACAAATATGA'+
      'AAAGTTACACTTTTCAGGTGT'+
      'TTTTG'+
      'AGAAATACAA'
    exp = exp.revcom

    graph.nodes[61].sequence.should == exp
  end
end
