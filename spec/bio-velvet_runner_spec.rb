require File.expand_path(File.dirname(__FILE__) + '/spec_helper')

describe "BioVelvet" do
  it "should run velvet" do
    begin
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
      exp.gsub!(/ /,'')

      File.open(contigs_path).read.should == exp

      File.exist?(result.contigs_path).should == true
      File.exist?(result.last_graph_path).should == true
      File.exist?(result.stats_path).should == true
    end
  end
end
