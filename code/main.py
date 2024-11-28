import re
from q_stat import *
from filter import *


def argparserLocal():
    from argparse import ArgumentParser
    '''Argument parser for the commands'''
    parser = ArgumentParser(prog='fastq', description='Filtering DNA reads by quality score, read length, and barcode. Then, create a new filtered FASTQ file.')

    subparsers = parser.add_subparsers(
        title='commands', description='Please choose command below:',
        dest='command'
    )
    subparsers.required = True

    #stat
    cal_q_command = subparsers.add_parser('calQ', help='Calculate mean of quality score for each read.')
    cal_q_command.add_argument("-f", "--file", type=str, default=None,
                            help="Provide FASTQ file or gzip file to program.")
    cal_q_command.add_argument("-b", "--barcode", action='store_true',
                            help="Calculate mean of quality score for each barcode.")
    cal_q_command.add_argument("-m", "--median", action='store_true',
                            help="Calculate median of quality score.")
    
    length_command = subparsers.add_parser('len', help=': Calculate length for each read.')
    length_command.add_argument("-f", "--file", type=str, default=None,
                            help="Provide FASTQ file or gzip file to program.")
    length_command.add_argument("-b", "--barcode", action='store_true',
                            help="Calculate mean of read length in each barcode.")
    length_command.add_argument("-m", "--median", action='store_true',
                            help="Calculate median of read length.")
    
    q_summary_command = subparsers.add_parser('Qsum', help='Calculate mean of quality score and length for each read.')
    q_summary_command.add_argument("-f", "--file", type=str, default=None,
                            help="Provide FASTQ file or gzip file to program.")
    q_summary_command.add_argument("-b", "--barcode", action='store_true',
                            help="Calculate mean of read length and mean of quality score in each barcode.")
    q_summary_command.add_argument("-m", "--median", action='store_true',
                            help="Calculate median of read length and median of quality score.")

    #filter
    filter_command = subparsers.add_parser('filter', help='Show filtering status for each read.')
    filter_command.add_argument("-f", "--file", type=str, default=None,
                            help="Provide FASTQ file or gzip file to program.")
    filter_command.add_argument("-l", "--length", type=int, default=100,
                            help="Filter reads by length, default length >=100 bp.")
    filter_command.add_argument("-q", "--qscore", type=int, default=30,
                            help="Filter reads by mean of quality score, default quality score >= 30.")
    filter_command.add_argument("-b", "--barcode", action='store_true',
                            help="Show filtering status for each barcode.")
    filter_command.add_argument("-m", "--median", action='store_true',
                            help="Calculate median of read length and median of quality score.")
    
    new_fastq_command = subparsers.add_parser('newFASTQ', help='Create a new FASTQ file with filtered reads.')
    new_fastq_command.add_argument("-f", "--file", type=str, default=None,
                            help="Provide FASTQ file or gzip file to program.")
    new_fastq_command.add_argument("-l", "--length", type=int, default=100,
                            help="Filter reads by length, default length >=100 bp.")
    new_fastq_command.add_argument("-q", "--qscore", type=int, default=30,
                            help="Filter reads by mean of quality score, default quality score >= 30.")
    new_fastq_command.add_argument("-m", "--median", action='store_true',
                            help="Calculate median of read length and median of quality score.")
    new_fastq_command.add_argument("-n", "--name", type=str, default= 'summary_data.fastq',
                            help="Name the new filtered file, default name = summary_data.fastq.")
    

    # parser.print_help()
    return parser

def read_zipfile(input_file):
    '''read data in compressed FASTQ file by gzip.'''    
    import gzip

    # เปิดไฟล์ .gz ด้วย read text mode (rt): อ่านข้อมูลในโหมดข้อความ (string) ไม่ใช่โหมดไบนารี bytes
    with gzip.open(input_file, "rt") as input_data:
        for line in input_data:
            yield line.strip()  # ส่งคืนบรรทัดโดยตัดช่องว่าง แบบ iterator

def read_file(input_file):
    '''read data in the FASTQ file (.fastq)'''

    with open(input_file, "r") as input_data:
        for line in input_data:
            yield line.strip()  # ส่งคืนบรรทัดโดยตัดช่องว่าง แบบ iterator

def get_data(input_file):
    '''get data from file'''
    if '.gz' in input_file :
        return read_zipfile(input_file)
    else:
        return read_file(input_file )

def line1(line):
    '''Search barcode and ID.'''
    return re.search(r"barcode=(\S*)",line).group(1), re.search(r"^@(\S*)",line).group(1)

def line4(line, command, option_median, filter_l = 100, filter_q = 30) :
    '''Select function from q_stat.py to work with line 4.'''
    if command == "calQ":
        if option_median == False:
            return calmean(line), SDquality(line)
        elif option_median == True:
            return calquantile(line, True)

    elif command == "len":
        return lenofreads(line)

    elif command == "Qsum":
        if option_median == False:
            return lenofreads(line), calmean(line), SDquality(line)
        elif option_median == True:
            return lenofreads(line), calquantile(line, True)

    elif command == "filter":
        if option_median == False:
            length = lenofreads(line)
            quality = calmean(line)
            return length, quality, status_read(length, quality, filter_l, filter_q)
        elif option_median == True:
            length = lenofreads(line)
            median = calquantile(line)
            return length, median, status_read(length, median, filter_l, filter_q)
        
    elif command == "newFASTQ":
        if option_median == False:
            return status_read(lenofreads(line), calmean(line), filter_l, filter_q)
        elif option_median == True:
            return status_read(lenofreads(line), calquantile(line), filter_l, filter_q)

def option_barcode(input_file, command, option_median, filter_l = 100, filter_q = 30):
    '''Collects data for each barcode.'''
    list_read = {}
    for i, line in enumerate(get_data(input_file)):
        if i % 4 == 0 :
            barcode = re.search(r"barcode=(\S*)",line).group(1)
            list_read.setdefault(barcode,{'num_of_read': 0,'q':[],'l':[],'status':[]})
        elif i % 4 == 3 :
            list_read[barcode]['num_of_read'] = list_read[barcode]['num_of_read'] + 1
            if command == "calQ":
                if option_median == False:
                    list_read[barcode]['q'].append(calmean(line)) 
                elif option_median == True:
                    list_read[barcode]['q'].append(calquantile(line)) 

            elif command == "len":
                    list_read[barcode]['l'].append(lenofreads(line))

            elif command == "Qsum":
                list_read[barcode]['l'].append(lenofreads(line))
                if option_median == False:
                    list_read[barcode]['q'].append(calmean(line))
                elif option_median == True:
                    list_read[barcode]['q'].append(calquantile(line))
            
            elif command == "filter":
                if option_median == False:
                    list_read[barcode]['status'].append(status_read(lenofreads(line), calmean(line), filter_l, filter_q)) 
                elif option_median == True:
                    list_read[barcode]['status'].append(status_read(lenofreads(line), calquantile(line), filter_l, filter_q))
    return list_read 

def status_read(length, quality, filter_l, filter_q):
    '''Summary status of each read.'''
    if filter_lenght(length, filter_l) == "Pass" and filter_quality(quality, filter_q)  == "Pass":
        return "Pass"
    else:
        return "Not pass"
    
def main():
    '''Main Function for run all feature'''
    parser = argparserLocal()
    args = parser.parse_args()

    input_file = args.file

    if args.command == 'calQ':
        if args.file == None:
            exit(parser.parse_args(['calQ','-h']))
        if args.barcode == False:
            if args.median == False:
                print(f"|                read id                 |   barcode   | Q score |   SD   |")
                print( "---------------------------------------------------------------------------")
                for i, line in enumerate(get_data(input_file)):
                    if i % 4 == 0 :
                        barcode, id = line1(line)
                    elif i % 4 == 3 :
                        quality, sd = line4(line, args.command, args.median)
                        print(f"|{id:^40s}|{barcode:^13s}|{quality:^9,.02f}|{sd:>6,.02f}  |")
            elif args.median == True:
                print(f"|                read id                 |   barcode   | median Q score |   IQR   |")
                print( "-----------------------------------------------------------------------------------")
                for i, line in enumerate(get_data(input_file)):
                    if i % 4 == 0 :
                        barcode, id = line1(line)
                    elif i % 4 == 3 :
                        median, IQR = line4(line, args.command, args.median)
                        print(f"|{id:^40s}|{barcode:^13s}|{median:^16,.02f}|{IQR:>7,.02f}  |")

        elif args.barcode == True:
            list_read = option_barcode(input_file, args.command, args.median)
            if args.median == False:
                print(f"|   barcode   | number of reads | mean Q score |   SD   |  min  |  max  |")
                print( "-------------------------------------------------------------------------")
                for barcode in list_read:
                    mean_q, sd, min_q, max_q = callenandbar(list_read[barcode],'q')
                    print(f"|{barcode:^13s}|{list_read[barcode]['num_of_read']:^17d}|{mean_q:^14,.02f}|{sd:>6,.02f}  |{min_q:^7,.02f}|{max_q:^7,.02f}|")
            elif args.median == True:
                print(f"|   barcode   | number of reads | median Q score |   IQR   |  min  |  max  |")
                print( "----------------------------------------------------------------------------")
                for barcode in list_read:
                    median, IQR, min_q, max_q = callenandbar_median(list_read[barcode],'q')
                    # print(list_read[barcode]['q'])
                    print(f"|{barcode:^13s}|{list_read[barcode]['num_of_read']:^17d}|{median:^16,.02f}|{IQR:>7,.02f}  |{min_q:^7,.02f}|{max_q:^7,.02f}|")

        
    if args.command == 'len':
        if args.file == None:
            exit(parser.parse_args(['len','-h']))
        if args.barcode == False:
            print(f"|                read id                 |   barcode   | length |")
            print( "-----------------------------------------------------------------")
            for i, line in enumerate(get_data(input_file)):
                if i % 4 == 0 :
                    barcode, id = line1(line)
                elif i % 4 == 1 :
                    length = line4(line, args.command, args.median)
                    print(f"|{id:^40s}|{barcode:^13s}|{length:>7d} |")
        elif args.barcode == True:
            list_read = option_barcode(input_file, args.command, args.median)
            if args.median == False:
                print(f"|   barcode   | number of reads | mean length |     SD     |   min   |   max   |")
                print( "--------------------------------------------------------------------------------")
                for barcode in list_read:
                    mean_l, sd_l, min_l, max_l = callenandbar(list_read[barcode],'l')
                    print(f"|{barcode:^13s}|{list_read[barcode]['num_of_read']:^17d}|{mean_l:>11,.02f}  |{sd_l:>10,.02f}  |{min_l:>7d}  |{max_l:>7d}  |")
            elif args.median == True:
                print(f"|   barcode   | number of reads | median length |     IQR     |   min   |   max   |")
                print( "-----------------------------------------------------------------------------------")
                for barcode in list_read:
                    median, IQR, min_l, max_l = callenandbar_median(list_read[barcode],'l')
                    print(f"|{barcode:^13s}|{list_read[barcode]['num_of_read']:^17d}|{median:>13,.02f}  |{IQR:>11,.02f}  |{min_l:>7d}  |{max_l:>7d}  |")
        

    if args.command == 'Qsum':
        if args.file == None:
            exit(parser.parse_args(['Qsum','-h']))
        if args.barcode == False:
            if args.median == False:
                print(f"|                read id                 |   barcode   | length | Q score | SD Q score |")
                print( "----------------------------------------------------------------------------------------")
                for i, line in enumerate(get_data(input_file)):
                    if i % 4 == 0 :
                        barcode, id = line1(line)
                    elif i % 4 == 3 :
                        length, quality, sd = line4(line, args.command, args.median)
                        print(f"|{id:^40s}|{barcode:^13s}|{length:>7d} |{quality:^9,.02f}|{sd:^12,.02f}|")
            elif args.median == True:
                print(f"|                read id                 |   barcode   | length | median Q score | IQR Q score |")
                print( "------------------------------------------------------------------------------------------------")
                for i, line in enumerate(get_data(input_file)):
                    if i % 4 == 0 :
                        barcode, id = line1(line)
                    elif i % 4 == 3 :
                        # print(line4(line, args.command, args.median))
                        length, (median, IQR) = line4(line, args.command, args.median)
                        print(f"|{id:^40s}|{barcode:^13s}|{length:>7d} |{median:^16,.02f}|{IQR:>9,.02f}    |")
        elif args.barcode == True:
            list_read = option_barcode(input_file, args.command, args.median)
            if args.median == False:
                print(f"|   barcode   | number of reads | mean length |  SD length  |  min length  |  max length  | mean Q score | SD Q score |  min Q score  |  max Q score  |")
                print( "-------------------------------------------------------------------------------------------------------------------------------------------------------")
                for barcode in list_read:
                    mean_q, sd_q, min_q, max_q = callenandbar(list_read[barcode],'q')
                    mean_l, sd_l, min_l, max_l = callenandbar(list_read[barcode],'l')  
                    print(f"|{barcode:^13s}|{list_read[barcode]['num_of_read']:^17d}|{mean_l:>11,.02f}  |{sd_l:>11,.02f}  |{min_l:>12d}  |{max_l:>12d}  |{mean_q:^14,.02f}|{sd_q:>8,.02f}    |{min_q:^15,.02f}|{max_q:^15,.02f}|")
            elif args.median == True:
                print(f"|   barcode   | number of reads | median length |  IQR length  |  min length  |  max length  | median Q score |  IQR Q score |  min Q score  |  max Q score  |")
                print( "--------------------------------------------------------------------------------------------------------------------------------------------------------------")
                for barcode in list_read:
                    median_q, IQR_q, min_q, max_q = callenandbar_median(list_read[barcode],'q')
                    median_l, IQR_l, min_l, max_l = callenandbar_median(list_read[barcode],'l')
                    print(f"|{barcode:^13s}|{list_read[barcode]['num_of_read']:^17d}|{median_l:>13,.02f}  |{IQR_l:>12,.02f}  |{min_l:>12d}  |{max_l:>12d}  |{median_q:^16,.02f}|{IQR_q:>10,.02f}    |{min_q:^15,.02f}|{max_q:^15,.02f}|")

    if args.command == 'filter':
        if args.file == None:
            exit(parser.parse_args(['filter','-h']))
        if args.barcode == False:
            if args.median == False:
                print(f"|                read id                 |   barcode   | length | Q score |  status  |")
                print( "--------------------------------------------------------------------------------------")
                for i, line in enumerate(get_data(input_file)):
                    if i % 4 == 0 :
                        barcode, id = line1(line)
                    elif i % 4 == 3 :
                        length, quality, status = line4(line, args.command, args.median, args.length, args.qscore)
                        print(f"|{id:^40s}|{barcode:^13s}|{length:>7d} |{quality:^9,.02f}|{status:^10s}|")
            elif args.median == True:
                print(f"|                read id                 |   barcode   | length | median Q score |  status  |")
                print( "---------------------------------------------------------------------------------------------")
                for i, line in enumerate(get_data(input_file)):
                    if i % 4 == 0 :
                        barcode, id = line1(line)
                    elif i % 4 == 3 :
                        length, median, status = line4(line, args.command, args.median, args.length, args.qscore)
                        print(f"|{id:^40s}|{barcode:^13s}|{length:>7d} |{median:^16,.02f}|{status:^10s}|")

        elif args.barcode == True:
            list_read = option_barcode(input_file, args.command, args.median, args.length, args.qscore)
            print(f"|   barcode   | number of reads | number of passed reads |  status  |")
            print( "---------------------------------------------------------------------")
            for barcode in list_read:
                read_pass, status_percent = filter_bar(list_read[barcode])
                status_percent = str(f"{status_percent:.02f}") + ' %'
                print(f"|{barcode:^13s}|{list_read[barcode]['num_of_read']:^17d}|{read_pass:^24d}|{status_percent:>9s} |")

                

    if args.command == 'newFASTQ':
        if args.file == None:
            exit(parser.parse_args(['newFASTQ','-h']))
        else: 
            if not re.search(r'\.f.*q',args.name):
                    args.name += '.fastq'
            file_name = name_file(args.name)
            list_read = {}
            print(f"|   barcode   | number of reads | number of passed reads |")
            print( "----------------------------------------------------------")
            for i, line in enumerate(get_data(input_file)):
                if i % 4 == 0 :
                    barcode = re.search(r"barcode=(\S*)",line).group(1)
                    list_read.setdefault(barcode,{'num_of_read': 0,'pass':0})
                    line_1 = line
                elif i % 4 == 1:
                    line_2 = line
                elif i % 4 == 3 :
                    line_4 = line
                    list_read[barcode]['num_of_read'] = list_read[barcode]['num_of_read'] + 1
                    status = line4(line, args.command, args.median, args.length, args.qscore)
                    if status == "Pass":
                        list_read[barcode]['pass'] = list_read[barcode]['pass'] + 1
                        write_file(line_1,line_2,line_4,file_name) 
            for barcode in list_read:
                print(f"|{barcode:^13s}|{list_read[barcode]['num_of_read']:^17d}|{list_read[barcode]['pass']:^24d}|")
            print( "----------------------------------------------------------\n")
            print("your new FASTQ file name:", file_name)
            print("\n")
           

if __name__ == '__main__':
    # main()
    input_file = 'data_test.fastq'
    # 'data_test_2.fastq'
    # 'data_test_3.fastq'
    # input_file = 'ont_reads.project.fastq.gz'
    # test(input_file)
    # get_data(input_file)
    
    