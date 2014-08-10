files = Array.new
files.push({ :filename => "SRR505743",    :reads => 225257463, :length => 101, :shift => 33 })
files.push({ :filename => "SRR505744",    :reads => 244674787, :length => 101, :shift => 33 })
files.push({ :filename => "SRR505745",    :reads => 226165495, :length => 101, :shift => 33 })
files.push({ :filename => "SRR505746",    :reads => 252345156, :length => 101, :shift => 33 })

files.push({ :filename => "SRR557711",    :reads => 36037705,  :length => 36,  :shift => 33 })
files.push({ :filename => "SRR557723",    :reads => 37525647,  :length => 36,  :shift => 33 })

files.push({ :filename => "SRR639080_1",  :reads => 16099716,  :length => 101, :shift => 33 })

files.push({ :filename => "SRR985867",    :reads => 21129136,  :length => 50,  :shift => 33 })

files.push({ :filename => "SRR988190",    :reads => 56746324,  :length => 202, :shift => 33 })
files.push({ :filename => "SRR988193",    :reads => 48330712,  :length => 202, :shift => 33 })

files.push({ :filename => "SRR1029924",   :reads => 87105048,  :length => 50,  :shift => 33 })
files.push({ :filename => "SRR1029925",   :reads => 83487348,  :length => 50,  :shift => 33 })

files.push({ :filename => "SRR1030717",   :reads => 87725913,  :length => 97,  :shift => 33 })
files.push({ :filename => "SRR1163160_1", :reads => 82517320,  :length => 100, :shift => 33 })
files.push({ :filename => "B1491_TAGCTT_L001_R1", :reads => 68679919, :length => 101, :shift => 33 })
files.push({ :filename => "B1492_AGTTCC_L001_R1", :reads => 97358268, :length => 101, :shift => 33 })


execPath = "/home/hedtke/git/SequenceTrimming/"
dataPath = "/space/GrosseHedtkeLemnianMuellerHannemann/"
timeCmd = "time -f 'total: %e \\t\\t user: %U' "
cacheCmd = "sync; echo 3 > /proc/sys/vm/drop_caches"
maxRepeats = 3
debug = false

## DISK SPEED
#for i in 1..maxRepeats
#    if !debug
#        system(cacheCmd)
#        print "!!! leere Cache !!!"
#        sleep 10
#    end
#    files.each do |file|
#        print "diskSpeed " + file[:filename].to_s + " run " + i.to_s + "\n"
#        mycmd = timeCmd.to_s + execPath.to_s + "tools_for_paper/diskSpeed "
#        mycmd += "-i " + dataPath.to_s + file[:filename].to_s + ".fastq "
#        mycmd += "-l " + file[:length].to_s + " "
#        mycmd += "-r " + file[:reads].to_s #+ " > /dev/null 2>&1'"
#        if debug
#            print mycmd + "\n"
#        else
#            system(mycmd)
#            print "\n\n\n"
#            sleep 5
#        end
#    end
#end
#
## 0-ZEROS
#for i in 1..maxRepeats
#    for t in [25,30]
#        if !debug
#            system(cacheCmd)
#            print "!!! leere Cache !!!"
#            sleep 10
#        end
#        files.each do |file|
#            print "0-zeros T=" + t.to_s + " " + file[:filename].to_s + " run " + i.to_s + "\n"
#            mycmd = timeCmd.to_s + execPath.to_s + "trimZeroOne "
#            mycmd += "-i " + dataPath.to_s + file[:filename].to_s + ".fastq "
#            mycmd += "-l " + file[:length].to_s + " "
#            mycmd += "-r " + file[:reads].to_s + " "
#            mycmd += "-s " + file[:shift].to_s + " "
#            mycmd += "-t " + t.to_s + " "
#            if debug
#                print mycmd + "\n"
#            else
#                system(mycmd)
#                print "\n\n\n"
#                sleep 5
#            end
#        end
#    end
#end
#
## 5-ZEROS
#for i in 1..maxRepeats
#    for t in [25,30]
#        if !debug
#            system(cacheCmd)
#            print "!!! leere Cache !!!"
#            sleep 10
#        end
#        files.each do |file|
#            print "5-zeros T=" + t.to_s + " " + file[:filename].to_s + " run " + i.to_s + "\n"
#            mycmd = timeCmd.to_s + execPath.to_s + "trimZeroOneZerosAllowed -z 5 "
#            mycmd += "-i " + dataPath.to_s + file[:filename].to_s + ".fastq "
#            mycmd += "-l " + file[:length].to_s + " "
#            mycmd += "-r " + file[:reads].to_s + " "
#            mycmd += "-s " + file[:shift].to_s + " "
#            mycmd += "-t " + t.to_s + " "
#            if debug
#                print mycmd + "\n"
#            else
#                system(mycmd)
#                print "\n\n\n"
#                sleep 5
#            end
#        end
#    end
#end

# 10-ZEROS
for i in 1..maxRepeats
    for t in [25,30]
        if !debug
            system(cacheCmd)
            print "!!! leere Cache !!!"
            sleep 10
        end
        files.each do |file|
            print "10-zeros T=" + t.to_s + " " + file[:filename].to_s + " run " + i.to_s + "\n"
            mycmd = timeCmd.to_s + execPath.to_s + "trimZeroOneZerosAllowed -z 10 "
            mycmd += "-i " + dataPath.to_s + file[:filename].to_s + ".fastq "
            mycmd += "-l " + file[:length].to_s + " "
            mycmd += "-r " + file[:reads].to_s + " "
            mycmd += "-s " + file[:shift].to_s + " "
            mycmd += "-t " + t.to_s + " "
            if debug
                print mycmd + "\n"
            else
                system(mycmd)
                print "\n\n\n"
                sleep 5
            end
        end
    end
end

## m-mean
#for i in 1..maxRepeats
#    for m in [25,30,35]
#        if !debug
#            system(cacheCmd)
#            print "!!! leere Cache !!!"
#            sleep 10
#        end
#        files.each do |file|
#            print "m-mean m=" + m.to_s + " " + file[:filename].to_s + " run " + i.to_s + "\n"
#            mycmd = timeCmd.to_s + execPath.to_s + "trimIntegerMean "
#            mycmd += "-i " + dataPath.to_s + file[:filename].to_s + ".fastq "
#            mycmd += "-l " + file[:length].to_s + " "
#            mycmd += "-r " + file[:reads].to_s + " "
#            mycmd += "-s " + file[:shift].to_s + " "
#            mycmd += "-m " + m.to_s + " "
#            if debug
#                print mycmd + "\n"
#            else
#                system(mycmd)
#                print "\n\n\n"
#                sleep 5
#            end
#        end
#    end
#end
#
## parallel 35-mean
#for i in 1..maxRepeats
#    for w in [2,4,8]
#        if !debug
#            system(cacheCmd)
#            print "!!! leere Cache !!!"
#            sleep 10
#        end
#        files.each do |file|
#            print "parallel 35-mean w=" + w.to_s + " " + file[:filename].to_s + " run " + i.to_s + "\n"
#            mycmd = timeCmd.to_s + execPath.to_s + "trimIntegerMean -m 35 "
#            mycmd += "-i " + dataPath.to_s + file[:filename].to_s + ".fastq "
#            mycmd += "-l " + file[:length].to_s + " "
#            mycmd += "-r " + file[:reads].to_s + " "
#            mycmd += "-s " + file[:shift].to_s + " "
#            mycmd += "-w " + w.to_s + " "
#            if debug
#                print mycmd + "\n"
#            else
#                system(mycmd)
#                print "\n\n\n"
#                sleep 5
#            end
#        end
#    end
#end
