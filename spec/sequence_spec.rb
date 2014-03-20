require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'

class String
  def revcom
    Bio::Sequence::NA.new(self).reverse_complement.to_s.upcase
  end
end

seq1_sequence = 'CACTTATCTCTACCAAAGATCACGATTTAGAATCAAACTATAAAGTTTTAGAAGATAAAG
TAACAACTTATACATGGGGATTCGGAGTTAAAAAAGTAGATTCAGAAAATATTTCAATAG
ATCTTGCAGGCGCAGCTTTTTCTGTTAGGGATAAAAATGGTAATGTAATTGGTAAATATA
CGTATGATTCTACTGGAAATGTGGTTTTATTAAAAGGAAAGGGTGTAACTGATAAAAATG
GACGAGTTATATTTACTGGTTTAAAAGAAGGAGATTACTTTATAAAAGAAGAAAAAGCTC
CTAAAGGGTATAGCCTTTTAAAAGAACCAGTAAAAGTTACTATAACAGCTCAAAAAGATG
ATAATGGAGAGTATACTGGTCAAGCAACTATATCTGTAACTAATGGCAATGAAGCTGGAA
GTATAATAAATAATATTACTATGAATGATGGCAATGTATTATTTAATGTACAAATTAAAA
ACTATGCTGGTATTTCACTTCCAGGTACAGG'.gsub("\n",'')
seq2_sequence = 'ACTGGTTTAAAAGAAGGAGATTACTTTATAAAAGAAGAAAAAGCTCCTAAAGGGTATAGC
CTTTTAAAAGAACCAGTAAAAGTTACTATAACAGCTCAAAAAGATGATAATGGAGAGTAT
ACTGGTCAAGCAACTATATCTGTAACTAATGGCAATGAAGCTGGAAGTATAATAAATAAT
ATTACTATGAATGATGGCAATGTATTATTTAATGTACAAATTAAAAACTATGCTGGTATT
TCACTTCCAGGTACAGGTGGAATTGGAACAGATGGATTCATTAAAATAGGGCTAGTTTTA
TTAGGGGTTGTTATTATTCTAGGTGCAGGATATGTTGTCTTAGATAAAAGAAAGAGAATT
TAATTAAATAAAAATATACTTCTTCTATTTTTAT'.gsub("\n",'')

describe "BioVelvet" do
  it "should be able to parse a whole Sequence file" do
    seqs = Bio::Velvet::Sequence.parse_from_file File.join(TEST_DATA_DIR, 'sequence_spec','Sequences')
    seqs.should be_kind_of(Bio::Velvet::Sequence)
    seqs.should be_kind_of(Hash)
    seqs.keys.should == [1,2]
    seqs[1].should == seq1_sequence
    seqs[2].should == seq2_sequence
  end

  it 'should be able to read in an interesting seq only' do
    seqs = Bio::Velvet::Sequence.parse_from_file File.join(TEST_DATA_DIR, 'sequence_spec','Sequences'),
    :interesting_read_ids => [1]
    seqs.keys.should == [1]
    seqs = Bio::Velvet::Sequence.parse_from_file File.join(TEST_DATA_DIR, 'sequence_spec','Sequences'),
    :interesting_read_ids => [2]
    seqs.keys.should == [2]
    seqs = Bio::Velvet::Sequence.parse_from_file File.join(TEST_DATA_DIR, 'sequence_spec','Sequences'),
    :interesting_read_ids => [1,2]
    seqs.keys.should == [1,2]
    seqs[1].should == seq1_sequence
  end

  it 'should be able to apply the grep hack' do
    seqs = Bio::Velvet::Sequence.parse_from_file File.join(TEST_DATA_DIR, 'sequence_spec','Sequences'),
    {:interesting_read_ids => [1], :apply_grep_hack => 500}
    seqs.keys.should == [1]
    seqs[1].should == seq1_sequence
  end

  it 'should be able to apply the grep hack when there is only just enough context' do
    seqs = Bio::Velvet::Sequence.parse_from_file File.join(TEST_DATA_DIR, 'sequence_spec','Sequences'),
    {:interesting_read_ids => [1], :apply_grep_hack => 9}
    seqs.keys.should == [1]
    seqs[1].should == seq1_sequence
  end

  it 'should warn when insufficient context is given' do
    expect {
      Bio::Velvet::Sequence.parse_from_file File.join(TEST_DATA_DIR, 'sequence_spec','Sequences'),
      {:interesting_read_ids => [1], :apply_grep_hack => 8}
      }.to raise_error
  end

  it 'should be able to handle multiple separated read ids with the grep hack' do
    s = Bio::Velvet::Sequence.parse_from_file File.join(TEST_DATA_DIR, 'sequence_spec','5seqs.fa.Sequences')
    s = ['']+s.values
    seq_file = File.join(TEST_DATA_DIR, 'sequence_spec','5seqs.fa.Sequences')

    seqs = Bio::Velvet::Sequence.parse_from_file seq_file,
    {:interesting_read_ids => [1], :apply_grep_hack => 9}
    seqs.keys.should == [1]
    seqs[1].should == s[1]

    seqs = Bio::Velvet::Sequence.parse_from_file seq_file,
    {:interesting_read_ids => [1,2], :apply_grep_hack => 9}
    seqs.keys.should == [1,2]
    seqs.values.should == s[1..2]

    seqs = Bio::Velvet::Sequence.parse_from_file seq_file,
    {:interesting_read_ids => [1,2], :apply_grep_hack => 2}
    seqs.keys.should == [1,2]
    seqs.values.should == s[1..2]

    seqs = Bio::Velvet::Sequence.parse_from_file seq_file,
    {:interesting_read_ids => [1,5], :apply_grep_hack => 2}
    seqs.keys.should == [1,5]
    seqs.values.should == [s[1],s[5]]

    seqs = Bio::Velvet::Sequence.parse_from_file seq_file,
    {:interesting_read_ids => [1,5], :apply_grep_hack => 3}
    seqs.keys.should == [1,5]
    seqs.values.should == [s[1],s[5]]

    seqs = Bio::Velvet::Sequence.parse_from_file seq_file,
    {:interesting_read_ids => [1,12], :apply_grep_hack => 3}
    seqs.keys.should == [1,12]
    seqs.values.should == [s[1],s[12]]

    expect {
      seqs = Bio::Velvet::Sequence.parse_from_file seq_file,
      {:interesting_read_ids => [12], :apply_grep_hack => 1}
      seqs.keys.should == [1,12]
      seqs.values.should == [s[1],s[12]]
      }.to raise_error
end
end
