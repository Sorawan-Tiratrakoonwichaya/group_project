from q_stat import calmean, lenofreads


def filter_quality(score,filter_scores=30):
    '''Filter by quality score'''
    if score >= filter_scores:
        return "Pass"
    else:
        return "Not pass"
    
def filter_lenght(length,filter_len=100):
    '''Filter by length'''
    if length >= filter_len:
        return "Pass"
    else:
        return "Not pass"
    
def filter_bar(dict):
    '''Calculate the percentage of passes per 1 barcode'''
    sumpass = 0
    total = 0
    for i in dict["status"]:
        total += 1
        if i == "Pass":
            sumpass += 1 
    # print(total, sumpass)
    return sumpass, (sumpass/total)*100

def name_file(new_file):
    '''Edit file name'''
    from pathlib import Path
    import re
    file_path = Path(new_file)
    if file_path.exists():
        new_file = file_path.stem + '1' + file_path.suffix
        while Path(new_file).exists():
            # print(re.search(r'([0-9]*)\.f',new_file).group(1), type(re.search(r'([0-9]*)\.f',new_file).group(1)))
            new_file = file_path.stem + str(int(re.search(r'([0-9]*)\.f',new_file).group(1))+1) + file_path.suffix
    return new_file


def write_file(Line1,Line2,Line4,new_file):
    '''write new file with filtered read'''
    with open(new_file,"a") as file:
        file.write(Line1+"\n")
        file.write(Line2+"\n")
        file.write("+\n")
        file.write(Line4+"\n")
    