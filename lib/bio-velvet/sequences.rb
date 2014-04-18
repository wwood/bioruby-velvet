#require 'hopcsv'
require 'bio'
require 'tempfile'

module Bio
  module Velvet
    # Parser and container class for textual Sequence files
    #
    # After parsing, the result is a hash of read_id => sequence
    # where read_id is an Integer and sequence a String
    #
    # The definition of this file is given in the velvet manual, at
    # http://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf
    class Sequences < Hash
      include Bio::Velvet::Logging

      def self.log
        self.new.log
      end

      # Options:
      # * :interesting_read_ids: If not nil, is a Set of nodes that we are interested in. Reads
      # not of interest will not be parsed in (the NR part of the velvet LastGraph file). Regardless all
      # nodes and edges are parsed in. Using this options saves both memory and CPU.
      # * :grep_hack: to make the parsing of read associations go even faster, a grep-based, rather
      # hacky method is applied to the graph file, so only sequence data of interesting_read_ids is presented
      # to the parser. This can save days of parsing time, but is a bit of a hack and its usage may
      # not be particularly future-proof. The value of this option is the amount of context coming out of grep
      # (the -A flag). In the Sequence file the sequences are wrapped at 60 characters, so you'll need at
      # least (longest_sequence_length / 60) + 2 amount of context. The reason for adding 2 is that the
      # parser will then be able to detect insufficient context and raise an Exception, without
      # throwing up false positive Exceptions.
      def self.parse_from_file(path_to_sequence_file, options={})
        seq_object = Bio::Velvet::Sequences.new

        if options[:apply_grep_hack]
          apply_grep_hack(seq_object, path_to_sequence_file, options[:interesting_read_ids], options[:apply_grep_hack])
        else
          # Parse all the sequences
          Bio::FlatFile.foreach(path_to_sequence_file) do |seq|
            read_id = seq.definition.split("\t")[1].to_i
            if options[:interesting_read_ids].nil? or options[:interesting_read_ids].include?(read_id)
              seq_object[read_id] = seq.seq.to_s
            end
          end
        end
        log.info "Read in #{seq_object.length} velvet stored sequences"
        return seq_object
      end

      private
      # Add the interesting sequences to the hash
      def self.apply_grep_hack(seq_object, path_to_sequence_file, interesting_read_ids, grep_context)
        return if interesting_read_ids.nil? or interesting_read_ids.empty?

        Tempfile.open('grep_v_hack') do |tempfile|
          # Create a file to pass to grep -f
          unless interesting_read_ids.nil?
            interesting_read_ids.each do |read_id|
              tempfile.puts "\t#{read_id}\t" #the read_id is the second field of the header
            end
          end
          tempfile.close

          cmd = "grep -F -A #{grep_context.inspect} -f #{tempfile.path} #{path_to_sequence_file.inspect}"
          # TODO: make this call more robust
          # grep_result = Bio::Commandeer.run cmd
          s, grep_result, stderr = systemu cmd

          # Parse the grepped out results
          current_read_id = nil
          current_seq = nil
          last_sequence_line_length = nil

          add_last_sequence = lambda do
            if current_read_id
              seq_object[current_read_id] = current_seq
            end
          end
          grep_result.each_line do |line|
            line.chomp!
            if line[0] == '>'
              # Process the last sequence
              add_last_sequence.call unless current_read_id.nil?

              # Assume the real sequence name contains no tabs
              read_id = line.split("\t")[1]
              raise "Unable to parse velvet Sequence file at this line #{line}" if read_id.nil?
              read_id = read_id.to_i
              if interesting_read_ids.include?(read_id)
                # if current_read_id is nil, then we know we are uninterested in this sequence
                current_read_id = read_id
              else
                current_read_id = nil
              end
              current_seq = nil
            elsif line == '--'
              # grep demarker.
              add_last_sequence.call unless current_read_id.nil?
              if last_sequence_line_length == 60
                raise "Parsing exception when parsing velvet Sequence file - grep hack too hacky. Sorry. Try modifying the code to increase the default amount of context grep is giving"
              end
            else
              # plain old sequence
              unless current_read_id.nil?
                current_seq ||= ''
                current_seq += line
                last_sequence_line_length = line.length
              end
            end
          end
          # process the last sequence
          if last_sequence_line_length == 60
            raise "Parsing exception when parsing velvet Sequence file at the end of the file - grep hack too hacky. Sorry. Try modifying the code to increase the default amount of context grep is giving"
          end
          add_last_sequence.call
        end
      end
    end
  end
end
