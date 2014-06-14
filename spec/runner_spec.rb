require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe "BioVelvet" do
  it "should run velvet" do
    reads_file = File.join TEST_DATA_DIR, 'runner_input.fa'
    result = Bio::Velvet::Runner.new.velvet(29,"-short #{reads_file}")
    result.should be_kind_of(Bio::Velvet::Result)

    contigs_path = File.join result.result_directory, 'contigs.fa'
    File.exist?(contigs_path).should == true
    exp = '>NODE_1_length_483_cov_1.474120
      CACTTATCTCTACCAAAGATCACGATTTAGAATCAAACTATAAAGTTTTAGAAGATAAAG
      TAACAACTTATACATGGGGATTCGGAGTTAAAAAAGTAGATTCAGAAAATATTTCAATAG
      ATCTTGCAGGCGCAGCTTTTTCTGTTAGGGATAAAAATGGTAATGTAATTGGTAAATATA
      CGTATGATTCTACTGGAAATGTGGTTTTATTAAAAGGAAAGGGTGTAACTGATAAAAATG
      GACGAGTTATATTTACTGGTTTAAAAGAAGGAGATTACTTTATAAAAGAAGAAAAAGCTC
      CTAAAGGGTATAGCCTTTTAAAAGAACCAGTAAAAGTTACTATAACAGCTCAAAAAGATG
      ATAATGGAGAGTATACTGGTCAAGCAACTATATCTGTAACTAATGGCAATGAAGCTGGAA
      GTATAATAAATAATATTACTATGAATGATGGCAATGTATTATTTAATGTACAAATTAAAA
      ACTATGCTGGTATTTCACTTCCAGGTACAGG'+"\n"
    exp.gsub!(/ / ,'')

    File.open(contigs_path).read.should == exp

    File.exist?(result.contigs_path).should == true
    File.exist?(result.last_graph_path).should == true
    File.exist?(result.stats_path).should == true
  end

  it 'should output to a real directory if told to do so' do
    reads_file = File.join TEST_DATA_DIR, 'runner_input.fa'
    outdir = '/tmp/my_velvet_assembly'
    File.exist?(outdir).should == false
    result = Bio::Velvet::Runner.new.velvet(29,"-short #{reads_file}",'',:output_assembly_path => outdir)
    result.should be_kind_of(Bio::Velvet::Result)
    result.result_directory.should == outdir
    File.exist?(result.result_directory).should == true
    File.exist?(outdir).should == true
    File.directory?(result.result_directory).should == true
    File.directory?(outdir).should == true

    contigs_path = File.join result.result_directory, 'contigs.fa'
    File.exist?(contigs_path).should == true
    exp = '>NODE_1_length_483_cov_1.474120
      CACTTATCTCTACCAAAGATCACGATTTAGAATCAAACTATAAAGTTTTAGAAGATAAAG
      TAACAACTTATACATGGGGATTCGGAGTTAAAAAAGTAGATTCAGAAAATATTTCAATAG
      ATCTTGCAGGCGCAGCTTTTTCTGTTAGGGATAAAAATGGTAATGTAATTGGTAAATATA
      CGTATGATTCTACTGGAAATGTGGTTTTATTAAAAGGAAAGGGTGTAACTGATAAAAATG
      GACGAGTTATATTTACTGGTTTAAAAGAAGGAGATTACTTTATAAAAGAAGAAAAAGCTC
      CTAAAGGGTATAGCCTTTTAAAAGAACCAGTAAAAGTTACTATAACAGCTCAAAAAGATG
      ATAATGGAGAGTATACTGGTCAAGCAACTATATCTGTAACTAATGGCAATGAAGCTGGAA
      GTATAATAAATAATATTACTATGAATGATGGCAATGTATTATTTAATGTACAAATTAAAA
      ACTATGCTGGTATTTCACTTCCAGGTACAGG'+"\n"
    exp.gsub!(/ / ,'')

    File.open(contigs_path).read.should == exp

    File.exist?(result.contigs_path).should == true
    File.exist?(result.last_graph_path).should == true
    File.exist?(result.stats_path).should == true

    #Clean up
    File.exist?(outdir).should == true
    FileUtils.remove_entry_secure(result.result_directory, true)
    File.exist?(outdir).should == false
  end

  it 'should detect binary_version' do
    Bio::Velvet::Runner.new.binary_version.match(/^1\..+/).should_not be_nil #will velvet every make it out of version 1.X ? If so, this test will fail.
  end
end
