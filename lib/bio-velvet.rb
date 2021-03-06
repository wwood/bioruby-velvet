require 'bio'
require 'bio-logger'
Bio::Log::LoggerPlus.new('bio-velvet')
module Bio::Velvet
  module Logging
    def log
      Bio::Log::LoggerPlus['bio-velvet']
    end
  end
end

require 'bio-velvet/graph'
require 'bio-velvet/runner'
require 'bio-velvet/sequences'
require 'bio-velvet/sequence_names'
