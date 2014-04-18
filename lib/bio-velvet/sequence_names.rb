require 'tempfile'
require 'bio-commandeer'

module Bio
  module Velvet
    # Methods for dealing with the CnyUnifiedSeq.names file output when
    # the -create_binary flag is set in velveth.
    class CnyUnifiedSeqNamesFile
      # Return a hash of seqname => (Array of CnyUnifiedSeqNamesFileEntry objects)
      # created by parsing the CnyUnifiedSeq.names file. Sometimes sequences
      # can be found multiple times e.g. if fwd and rev of a pair is delineated
      # after a space in the input names.
      def self.extract_entries(path_to_cny_unified_seq_names_file, entry_names)
        # Create results hash
        to_return = {}
        entry_names.each do |name|
          to_return[name] = []
        end

        Hopcsv.foreach(path_to_cny_unified_seq_names_file,"\t") do |row|
          name = row[0][1...row.length] #remove '>' at the start of the name
          next unless to_return.key?(name) #ignore uninsteresting sequences

          entry = CnyUnifiedSeqNamesFileEntry.new
          entry.name = name
          entry.read_id = row[1].to_i
          entry.category = row[2].to_i
          to_return[name].push entry
        end
        return to_return
      end

      # These files can be quite big, so this method
      def self.extract_entries_using_grep_hack(path_to_cny_unified_seq_names_file, entry_names)
        to_return = nil
        Tempfile.open('velvet_names_grep_hack_in') do |input|
          entry_names.each do |name|
            input.puts ">#{name}\t"
          end
          input.close #flush

          Tempfile.open('velvet_names_grep_hack_result') do |output|
            Bio::commandeer.run "grep -F -f #{input.path} #{path_to_cny_unified_seq_names_file.inspect} >#{output.path}"
            to_return = extract_entries output.path, entry_names
          end
        end
      end
    end

    class CnyUnifiedSeqNamesFileEntry
      attr_accessor :name, :read_id, :category
    end
  end
end