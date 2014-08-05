files = Array.new
files.push({ :filename => "SRR505743", :reads => 225257463, :length => 101, :shift => 33 })
files.push({ :filename => "SRR505744", :reads => 244674787, :length => 101, :shift => 33 })
files.push({ :filename => "SRR505745", :reads => 226165495, :length => 101, :shift => 33 })
files.push({ :filename => "SRR505746", :reads => 252345156, :length => 101, :shift => 33 })
files.push({ :filename => "SRR985867", :reads => 21129136, :length => 50, :shift => 33 })
files.push({ :filename => "SRR988190", :reads => 56746324, :length => 202, :shift => 33 })
files.push({ :filename => "SRR988193", :reads => 48330712, :length => 202, :shift => 33 })
files.push({ :filename => "SRR1029924", :reads => 87105048, :length => 50, :shift => 33 })
files.push({ :filename => "SRR1029925", :reads => 83487348, :length => 50, :shift => 33 })
files.push({ :filename => "SRR1030717", :reads => 87725913, :length => 97, :shift => 33 })

for i in 1..5
    print 'sh -c "sync; echo 3 > /proc/sys/vm/drop_caches"'
    #system(sh -c "sync; echo 3 > /proc/sys/vm/drop_caches")
    files.each do |file|
        mycmd = "time /home/hedtke/SequenceTrimming/tools_for_paper/diskSpeed "
        mycmd += "-i /data/hedtke/" + file[:filename].to_s + ".fastq "
        mycmd += "-l " + file[:length].to_s + " "
        mycmd += "-r " + file[:reads].to_s + " "
        mycmd += ">> /data/hedtke/diskSpeed_" + file[:filename].to_s + "_run_" + i.to_s + ".txt"
        print mycmd + "\n"
    end
end
