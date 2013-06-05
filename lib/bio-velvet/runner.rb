require 'tmpdir'
require 'systemu'

module Bio
  module Velvet
    class Runner
      # Run velveth and then velvetg, with the given kmer size
      #
      # The velveth_options hash must contain one of these keys:
      # 1. :short
      # 2. :short2
      # 3. :long
      # The value must a Array of paths to input sequence files.
      #
      # For binary flags, use true as the value in the hash
      def velvet(kmer_length, velveth_options, velvetg_options={})
        velveth_options_string = ''
        velveth_options.each do |key, value|
          if value == true
            velveth_options_string += ' '+key.to_s
          elsif key.to_s.match(/^short/) or key.to_s == 'long' or key.to_s == 'longPaired' or key.to_s == 'reference'
            # Multiple values allowed here
            velveth_options_string += ' ' + key.to_s
            if value.kind_of?(Array)
              velveth_options_string += value.join(' ')
            else
              velveth_options_string += value
            end
          else
            velveth_options_string += key.to_s+' '+ value
          end
        end
        velvetg_options_string = ''
        velvetg_options.each do |key, value|
          if value == true
            velvetg_options_string += ' '+key.to_s
          else
            velvetg_options_string += key.to_s+' '+ value
          end
        end

        # Run velveth
        status, stdout, stderr = systemu 'velveth'+opts
        stderr.should eq("")
    status.exitstatus.should eq(0)
    stdout.should eq("hello world")
      end
    end

    class VelvetRunnerException < Exception; end

    class Result
      attr_accessor :result_directory

      attr_a
    end
  end
end
