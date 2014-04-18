require File.expand_path(File.dirname(__FILE__) + '/spec_helper')
require 'bio'

include Bio::Velvet

describe "SeqeunceNames" do
  it 'should parse a whole file' do
    string = <<EOF
>read1	1	0
>read2	2	0
EOF
    Tempfile.open('test') do |tempfile|
      tempfile.print string
      tempfile.close

      names = Bio::Velvet::CnyUnifiedSeqNamesFile.extract_entries(
        tempfile.path,
        %w(read1 read2)
        )
      names.keys.should == %w(read1 read2)
      names['read1'].kind_of?(Array).should == true
      names['read1'].length.should == 1
      names['read1'][0].kind_of?(Bio::Velvet::CnyUnifiedSeqNamesFileEntry).should == true
      names['read1'].collect{|e| e.name}.should == ['read1']
      names['read1'].collect{|e| e.read_id}.should == [1]
      names['read1'].collect{|e| e.category}.should == [0]
      names['read2'].collect{|e| e.name}.should == ['read2']
      names['read2'].collect{|e| e.read_id}.should == [2]
      names['read2'].collect{|e| e.category}.should == [0]
    end
  end

  it 'should handle the grep hack' do
    string = <<EOF
>read1	1	0
>read2	2	0
EOF
    Tempfile.open('test') do |tempfile|
      tempfile.print string
      tempfile.close

      names = Bio::Velvet::CnyUnifiedSeqNamesFile.extract_entries_using_grep_hack(
        tempfile.path,
        %w(read1 read2)
        )
      names.keys.should == %w(read1 read2)
      names['read1'].kind_of?(Array).should == true
      names['read1'].length.should == 1
      names['read1'][0].kind_of?(Bio::Velvet::CnyUnifiedSeqNamesFileEntry).should == true
      names['read1'].collect{|e| e.name}.should == ['read1']
      names['read1'].collect{|e| e.read_id}.should == [1]
      names['read1'].collect{|e| e.category}.should == [0]
      names['read2'].collect{|e| e.name}.should == ['read2']
      names['read2'].collect{|e| e.read_id}.should == [2]
      names['read2'].collect{|e| e.category}.should == [0]

      names = Bio::Velvet::CnyUnifiedSeqNamesFile.extract_entries_using_grep_hack(
        tempfile.path,
        %w(read2)
        )
      names.keys.should == %w(read2)
      names['read2'].kind_of?(Array).should == true
      names['read2'].length.should == 1
      names['read2'][0].kind_of?(Bio::Velvet::CnyUnifiedSeqNamesFileEntry).should == true
      names['read2'].collect{|e| e.name}.should == ['read2']
      names['read2'].collect{|e| e.read_id}.should == [2]
      names['read2'].collect{|e| e.category}.should == [0]
    end
  end
end
