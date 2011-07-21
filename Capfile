require 'catpaws'

#generic settings
set :aws_access_key,  ENV['AMAZON_ACCESS_KEY']
set :aws_secret_access_key , ENV['AMAZON_SECRET_ACCESS_KEY']
set :ec2_url, ENV['EC2_URL']
set :ssh_options, { :user => "ubuntu", :keys=>[ENV['EC2_KEYFILE']]}
set :key, ENV['EC2_KEY']
set :key_file, ENV['EC2_KEYFILE']
set :ami, `curl http://mng.iop.kcl.ac.uk/cass_data/buckley_ami/AMIID`.chomp
set :instance_type,  'c1.xlarge'
set :working_dir, '/mnt/work'
set :nhosts, 1
set :group_name, '<PROJECT_NAME>'
set :snap_id, `cat SNAPID`.chomp 
set :vol_id, `cat VOLUMEID`.chomp 
set :ebs_size, 60 
set :availability_zone, 'eu-west-1a'
set :dev, '/dev/sdf'
set :mount_point, '/mnt/data'


# Allow Capfile.local to override these settings
begin
 load("Capfile.local")
rescue Exception
end

#start EC2 instances
#cap EC2:start

#make a new EBS volume from this snap 
#cap EBS:create

#attach your EBS
#cap EBS:attach

#format a new EBS
#cap EBS:format_xfs

#mount the EBS
#cap EBS:mount_xfs


#### Data Upload and Unpacking

desc "Upload data files"
task :upload_data, :roles => group_name do
  # TO DO
end
before 'upload_data', 'EC2:start'


# remove spaces in filenames - you'll need to run this before
# and after unzipping the files.
desc "Remove spaces in filenames"
task :rename_files, :roles => group_name do
  files = capture("ls #{mount_point}/*").split("\n")
  files.each do |f|
     newname = f.gsub(/\s+/, '_')
     newname = newname.gsub(/[\(\)\[\]\{\}\<\>]/,'_')
     if newname != f
       run "mv '#{f}' '#{newname}'"
     end
  end
end
before 'rename_files', 'EC2:start'


# zip hits probs with big files, use p7zip instead
desc "unzip data files"
task :unzip_files, :roles => group_name do
  run "sudo apt-get install -y unzip p7zip-full"
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/\.zip$/)}
  
  files.each do |f|
    run "cd #{mount_point} && 7za e '#{f}'"
  end
end
before "unzip_files", "EC2:start"





#### Quality Control

# convert the export.txt file to fastq for alignment
task :make_fastq, :roles => group_name do 
  upload("scripts/export2fastq.pl", "#{working_dir}/export2fastq.pl")
  run "chmod +x #{working_dir}/export2fastq.pl"
  run "sudo mv #{working_dir}/export2fastq.pl /usr/local/bin/"

  files = capture("ls #{mount_point}/*export.txt").split("\n")
  files = files.map {|f| f.chomp}

  files.each{|infile| 
    outfile = infile.sub('.txt', '.fastq')
    run "export2fastq.pl #{infile} > #{outfile}"
  } 

end
before 'make_fastq', 'EC2:start' 

# Run fastQC on the data files
desc "run fastqc"
task :fastqc, :roles => group_name do
  files = capture("ls #{mount_point}/*.fastq").split("\n")
  files = files.map {|f| f.chomp}
   
  files.each{|infile| 
    run "fastqc --outdir #{mount_point} #{infile}"
  } 

end
before 'fastqc', 'EC2:start'

# Pull the results back to the mng.iop.kcl.ac.uk server
desc "download fastqc files"
task :get_fastqc, :roles => group_name do
  `rm -Rf results/fastqc` #remove previous results
  `mkdir -p results/fastqc`
  files = capture "ls #{mount_point}/*fastqc.zip"
  files = files.split("\n")
  files.each{|f|
    outfile = f.sub(/.*\//,'')
    download( "#{f}", "results/fastqc/#{outfile}")
    `cd results/fastqc && unzip #{outfile} && rm #{outfile}`
  }
end
before "get_fastqc", 'EC2:start'



#### Alignment


#get the current mouse genome (which I already have on S3).
task :fetch_genome, :roles => group_name do
  run "mkdir -p #{working_dir}/indexes"
  run "cd #{working_dir}/indexes && curl https://s3-eu-west-1.amazonaws.com/genome-mm9-bowtie/mm9.ebwt.zip > mm9.ebwt.zip"
  run "rm -Rf #{working_dir}/indexes/chr*"
  run "cd  #{working_dir}/indexes && unzip -o mm9.ebwt.zip"
#?  run "export BOWTIE_INDEXES='#{working_dir}/indexes'"
end
before "fetch_genome","EC2:start"



# run bowtie on the fastq file
# This is recent illumina data, quals should be post v1.3
task :run_bowtie, :roles => group_name do

  files = capture("ls #{mount_point}/*.fastq").split("\n")
  files = files.map {|f| f.chomp}

  files.each{|infile|
    outfile = infile.sub('.fastq', '.sam')
    run("BOWTIE_INDEXES='#{working_dir}/indexes' && bowtie  --sam --best -k1 -l15 -n1 -m3 -p20 --solexa1.3-quals --chunkmbs 256  -q mm9 --quiet  #{infile}  > #{outfile} ")
  } 

end
before "run_bowtie", "EC2:start"


# Make binary BAM files from SAM
desc "make bam from sam"
task :to_bam, :roles => group_name do
  run "curl 'http://github.com/cassj/my_bioinfo_scripts/raw/master/genomes/mm9_lengths' > #{working_dir}/mm9_lengths"
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/\.sam$/)}
  files.each{|f| 
    f_out = f.sub('.sam', '.bam')
    run "samtools view -bt #{working_dir}/mm9_lengths -o #{mount_point}/#{f_out} #{mount_point}/#{f}"
  }
end
before "to_bam", "EC2:start"

# Sort the BAM files
desc "sort bam"
task :sort_bam, :roles => group_name do
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/\.bam/)}
  files.each{|f| 
    f_out = f.sub('.bam', '_sorted')
    run "samtools sort #{mount_point}/#{f}  #{mount_point}/#{f_out}"
  }
end
before "sort_bam", "EC2:start"


# Remove PCR Duplicate Reads
desc "remove duplicates"
task :rmdups, :roles => group_name do
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted\.bam/)}
  files.each{|f| 
    f_out = f.sub('sorted', 'sorted_nodups')
    run "cd #{mount_point} && samtools rmdup -s #{f}  #{f_out}"
  }
end
before "rmdups", "EC2:start"



# Index the BAM files
desc "index bam files"
task :index, :roles => group_name do
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted_nodups\.bam/)}
  files.each{|f| 
    f_out = f.sub('.bam', '.bai')
    run "cd #{mount_point} && samtools index  #{f} #{f_out}"
  }
end
before "index", "EC2:start"

# Create a summary of the files
desc "create a summary of the bam files"
task :flagstat, :roles => group_name do
 files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted_nodups\.bam/)}
  files.each{|f|
    f_out = f.sub('.bam', '.summary')
    run "cd #{mount_point} && samtools flagstat #{f} > #{f.out}"
  }

end
before "flagstat", "EC2:start"


# Pull the BAM files back to the mng.iop.kcl.ac.uk server
desc "download bam files"
task :get_bam, :roles => group_name do
  `rm -Rf results/alignment/bowtie` #remove previous results
  `mkdir -p results/alignment/bowtie`
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted_nodups/)}
  files.each{|f|
    download( "#{mount_point}/#{f}", "results/alignment/bowtie/#{f}")
  }
end
before "get_bam", "EC2:start"


# Need to do this or macs will run out of space to make WIG files on a 30GB partition.
# If you want the intermediate files, use a bigger EBS vol
desc "cleanup intermediate files"
task :bam_tidy, :roles =>group_name do
  files = capture "ls #{mount_point}"
  
  #don't delete the macs results, original export.txt files or final bam,bai files
  files = files.split("\n")
  files = files.reject{ |f| 
            f.match(/export.bam$/)   ||
            f.match(/^macs/)         ||  
            f.match(/sorted_nodups/) || 
            f.match(/export.txt$/)   ||
            f.match(/.*fastq.*/)     ||
            f.match(/README/)

          }
  files.each{|f| run("cd #{mount_point} && rm -Rf #{f}")}
end 
before "bam_tidy", 'EC2:start'





######## Peak Finding

#Need to be able to specify multiple treatment and control pairs here.

task :run_macs, :roles => group_name do

  treatment = "#{mount_point}/<TREATMENT FILE>"
  control = "#{mount_point}/<CONTROL FILE>"
  genome = 'mm'
  bws = [300]
  pvalues = [0.00001]

  #unsure what p values and bandwidths are appropriate, try a few?
  bws.each {|bw|
    pvalues.each { |pvalue|

      dir = "#{mount_point}/macs_#{bw}_#{pvalue}"
      run "rm -Rf #{dir}"
      run "mkdir #{dir}"

      macs_cmd =  "macs --treatment #{treatment} --control #{control} --name #{group_name} --format BAM --gsize #{genome} --bw #{bw} --pvalue #{pvalue}"
      run "cd #{dir} && #{macs_cmd}"
      
      dir = "#{mount_point}/macs_#{bw}_#{pvalue}_subpeaks"
      run "rm -Rf #{dir}"
      run "mkdir #{dir}"

      # With SubPeak finding
      # this will take a lot longer as you have to save the wig file
      macs_cmd =  "macs --treatment #{treatment} --control #{control} --name #{group_name} --format BAM --gsize #{genome} --call-subpeaks  --wig --bw #{bw} --pvalue #{pvalue}"
      run "cd #{dir} && #{macs_cmd}"

    }
  }
  
end
before 'run_macs', 'EC2:start'




# SICER

desc "bamToBed"
task :bamToBed, :roles => group_name do
   files = capture("ls #{mount_point}/*sorted_nodups.bam").split("\n")
   files = files.map {|f| f.chomp}
   files.each{|infile|
     f_out = infile.sub('.bam', '.bed')
     run "bamToBed -i #{infile} > #{f_out}"
   }

end
before 'bamToBed', 'EC2:start'


# this will need rewritten for your files.
desc "run SICER"
task :run_SICER, :roles => group_name do
   input = 'C18_input_CME142_s_7_export_sorted_nodups.bed'
   ip    = 'CME140_s_5_export_sorted_nodups.bed'
   

   species = 'mm9'
   thresh = 1
   window_size = 200
   fragment_size = 300
   effective_genome_fraction = '0.75' 
   gap_size = 600
   FDR = '0.1'

# /usr/local/bin/SICER [InputDir] [bed file] [control file] [OutputDir] [Species] [redundancy threshold] [window size (bp)] [fragment size] [effective genome fraction] [gap size (bp)] [FDR]

   run "mkdir -p #{mount_point}/SICER"
   run "SICER #{mount_point} #{ip} #{input} #{mount_point}/SICER #{species} #{thresh} #{window_size} #{fragment_size} #{effective_genome_fraction} #{gap_size} #{FDR}"






end
before 'run_SICER', 'EC2:start'


desc "convert peaks to IRanges"
task :peaks_to_iranges, :roles => group_name do
  run "cd #{working_dir} && rm -f peaksBed2IRanges.R"
  upload('scripts/peaksBed2IRanges.R' , "#{working_dir}/peaksBed2IRanges.R")
  run "cd #{working_dir} && chmod +x peaksBed2IRanges.R"

  macs_dirs = capture "ls #{mount_point}"
  macs_dirs = macs_dirs.split("\n")
  macs_dirs = macs_dirs.select { |d| d =~ /.*macs.*/ }
  macs_dirs.each{|d|
    d.chomp
    unless d.match('tgz')
      xlsfiles = capture "ls #{mount_point}/#{d}/*peaks.xls"
      xlsfiles = xlsfiles.split("\n")
      xlsfiles.each{|x|
         run "cd #{mount_point}/#{d} && Rscript #{working_dir}/peaksBed2IRanges.R #{x}"
      }
    end
  }
end
before 'peaks_to_iranges', 'EC2:start'

desc "annotate IRanges"
task :annotate_peaks, :roles => group_name do
  run "cd #{working_dir} && rm -f mm9RDtoGenes.R"
  upload('scripts/mm9RDtoGenes.R', "#{working_dir}/mm9RDtoGenes.R")
  run "cd #{working_dir} && chmod +x mm9RDtoGenes.R"
  macs_dirs = capture "ls #{mount_point}"
  macs_dirs = macs_dirs.split("\n").select { |d| d =~ /.*macs.*/ }
  macs_dirs.each{|d|
    unless d.match('tgz')
      rdfiles = capture "ls #{mount_point}/#{d}/*RangedData.RData"
      rdfiles = rdfiles.split("\n")
      rdfiles.each{|rd|
        unless rd.match(/negative/)
          run "cd #{working_dir} && Rscript #{working_dir}/mm9RDtoGenes.R #{rd}"
        end
      }
    end
  }
         
end
before 'annotate_peaks', 'EC2:start'

task :pack_peaks, :roles => group_name do
  macs_dirs = capture "ls #{mount_point}"
  macs_dirs = macs_dirs.split("\n").select {|f| f.match(/.*macs.*/)}
  macs_dirs.each{|d|
    unless d.match('tgz')
      run "cd #{mount_point} &&  tar --exclude *_wiggle* -cvzf #{d}.tgz #{d}"
    end
  }
  
end
before 'pack_peaks','EC2:start' 

task :get_peaks, :roles => group_name do
  macs_files = capture "ls #{mount_point}"
  macs_files = macs_files.split("\n").select {|f| f.match(/.*macs.*\.tgz/)}
  res_dir = 'results/alignment/bowtie/peakfinding/macs'
  `rm -Rf #{res_dir}`
  `mkdir -p #{res_dir}`
  macs_files.each{|f| 
    download("#{mount_point}/#{f}", "#{res_dir}/#{f}") 
    `cd #{res_dir} && tar -xvzf #{f}`
  }

end
before 'get_peaks', 'EC2:start'







#if you want to keep the results
#cap EBS:snapshot


#and then shut everything down:

# cap EBS:unmount
# cap EBS:detach
# cap EBS:delete - unless you're planning to use it again.
# cap EC2:stop

   




 
