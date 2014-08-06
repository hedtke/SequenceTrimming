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
    #system(sync; echo 3 > /proc/sys/vm/drop_caches)
    files.each do |file|
        print "diskSpeed " + file[:filename].to_s + " run " + i.to_s + "\n"
        mycmd = "time -f 'total: %e \t\t user: %U' /home/hedtke/git/SequenceTrimming/tools_for_paper/diskSpeed "
        mycmd += "-i /space/GrosseHedtkeLemnianMuellerHannemann/" + file[:filename].to_s + ".fastq "
        mycmd += "-l " + file[:length].to_s + " "
        mycmd += "-r " + file[:reads].to_s #+ " > /dev/null 2>&1'"
        #print mycmd + "\n"
        system(mycmd)
        print "\n\n\n"
        sleep 5
    end
end
