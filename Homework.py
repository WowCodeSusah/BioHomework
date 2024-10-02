# Function for Question 1
def Translate_DNA(dna_seq):
    # Makes all string case to uppercase to reduce errors
    seq = dna_seq.upper()

    # Checks if dna sequence is viable for codon search by dividing by 3
    checker = True
    if (len(dna_seq) % 3) == 0:
        checker = True
    else:
        checker = False

    # Checks if dna sequence is valid based on string characters
    valid = dna_seq.count("A") + dna_seq.count("G") + dna_seq.count("T") + dna_seq.count("C")
    if len(dna_seq) == valid:
        checker = True
    else:
        checker = False

    # Creates the complement by matching all the characters with their opposite character
    Complement = ""
    for x in seq:
        if x == 'A':
            Complement = Complement + "T"
        elif x == 'T':
            Complement = Complement + "A"
        elif x == 'C':
            Complement = Complement + "G"
        elif x == 'G':
            Complement = Complement + "C"
    
    # Using string replacement to replace all T's with U to find mRNA
    mRna = Complement.replace("T", "U")

    # Acid List to find keep all the codons
    acid = []

    # Codon Table without small letters
    codon_table = {
        'AUG': ('Methionine', 'M'), 'UUU': ('Phenylalanine', 'F'), 'UUC': ('Phenylalanine', 'F'),'UUA': ('Leucine', 'L'), 'UUG': ('Leucine', 'L'), 'UCU': ('Serine', 'S'), 'UCC': ('Serine', 'S'),
        'UCA': ('Serine', 'S'), 'UCG': ('Serine', 'S'), 'UAU': ('Tyrosine', 'Y'), 'UAC': ('Tyrosine', 'Y'),'UGU': ('Cysteine', 'C'), 'UGC': ('Cysteine', 'C'), 'UGG': ('Tryptophan', 'W'), 'CUU': ('Leucine', 'L'),
        'CUC': ('Leucine', 'L'), 'CUA': ('Leucine', 'L'), 'CUG': ('Leucine', 'L'), 'CCU': ('Proline', 'P'),'CCC': ('Proline', 'P'), 'CCA': ('Proline', 'P'), 'CCG': ('Proline', 'P'), 'CAU': ('Histidine', 'H'),'CAC': ('Histidine', 'H'), 'CAA': ('Glutamine', 'Q'), 'CAG': ('Glutamine', 'Q'), 'CGU': ('Arginine', 'R'),
        'CGC': ('Arginine', 'R'), 'CGA': ('Arginine', 'R'), 'CGG': ('Arginine', 'R'), 'AUU': ('Isoleucine', 'I'),'AUC': ('Isoleucine', 'I'), 'AUA': ('Isoleucine', 'I'), 'GUU': ('Valine', 'V'), 'GUC': ('Valine', 'V'),
        'GUA': ('Valine', 'V'), 'GUG': ('Valine', 'V'), 'GCU': ('Alanine', 'A'), 'GCC': ('Alanine', 'A'),'GCA': ('Alanine', 'A'), 'GCG': ('Alanine', 'A'), 'GAU': ('Aspartic acid', 'D'), 'GAC': ('Aspartic acid', 'D'),
        'GAA': ('Glutamic acid', 'E'), 'GAG': ('Glutamic acid', 'E'), 'GGU': ('Glycine', 'G'), 'GGC': ('Glycine', 'G'),'GGA': ('Glycine', 'G'), 'GGG': ('Glycine', 'G'), 'UAA': ('Stop', None), 'UAG': ('Stop', None),
        'UGA': ('Stop', None), 'ACU': ('Threonine', 'T'), 'ACC': ('Threonine', 'T'), 'ACA': ('Threonine', 'T'),'ACG': ('Threonine', 'T'), 'AAU': ('Asparagine', 'N'), 'AAC': ('Asparagine', 'N'), 'AAA': ('Lysine', 'K'),
        'AAG': ('Lysine', 'K'), 'AGU': ('Serine', 'S'), 'AGC': ('Serine', 'S'), 'AGA': ('Arginine', 'R'),'AGG': ('Arginine', 'R')
    }

    # Protein Letter table to show the letters in the amino acid
    protein_letters = {
        "A": "Ala","C": "Cys","D": "Asp","E": "Glu","F": "Phe","G": "Gly","H": "His","I": "Ile","K": "Lys","L": "Leu","M": 
        "Met","N": "Asn","P": "Pro","Q": "Gln","R": "Arg","S": "Ser","T": "Thr","V": "Val","W": "Trp","Y": "Tyr",
    }

    # Function to seperate all the codon bases into three parts
    fullAcid = ""
    acidCount = 0
    for currentAcid in mRna:
        fullAcid = fullAcid + currentAcid
        acidCount = acidCount + 1
        if acidCount == 3:
            acid.append(fullAcid)
            acidCount = 0
            fullAcid = ""

    # Function to search for all codon bases to find all the coresponding amino acid and letters
    resultText = "Aminoacid = "
    for acidCodon in acid:
        for codonSearch in codon_table:
            if acidCodon == codonSearch:
                for codonLetter in protein_letters:
                    if codon_table[codonSearch][1] == codonLetter:
                        resultText = resultText + " " + protein_letters[codonLetter] + " (" + codon_table[codonSearch][1] + ")"
        resultText = resultText + " -"
    resultText = resultText[:-1]

    # Prints the comploment, mRNA, and final result
    if checker == True:
        print ("Input DNA = " + seq)
        print("")
        print("Complement = " +Complement)
        print("mRNA = " + mRna)
        print(resultText)

    # Prints error message if any checks fail
    if checker == False:
        print("DNA sequence is Invalid")

# Function for Question 2
def frequencyCodon(aminoAcid):
    # Protein Letter table to show the letters in the amino acid
    protein_letters = {
        "A": ["GCC", "GCA", "GCG", "GCU"],"C": ["UGU", "UGC"],"D": ["GAU", "GAC"],"E": ["GAA", "GAG"],
        "F": ["UUU", "UUC"],"G": ["GGU", "GGC", "GGA", "GGG"],"H": ["CAU", "CAC"],"I": ["AU", "AUC", "AUA"],
        "K": ["AAA", "AAG"],"L": ["CUU", "CUC", "CUA", "CUG", "UUA", "UUG"],"M": ["AUG"],"N": ["AAU", "AAC"],
        "P": ["CCU", "CCC", "CCA", "CCG"],"Q": ["CAA", "CAG"],"R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],"T": ["ACU", "ACC", "ACA", "ACG"],"V": ["GUU", "GUC", "GUA", "GUG"],
        "W": ["UGG"],"Y": ["UAU", "UAC"],
    }

    # Double Checks for Valid amino acids
    for letters in aminoAcid.upper():
        if letters not in protein_letters:
            return print("Amino Acid is Invalid")

    TotalCombination = 1
    # checks all total amount of possible combinations
    for totalCom in aminoAcid.upper():
        for letter in protein_letters:
            if totalCom == letter:
                TotalCombination = TotalCombination * len(protein_letters[letter])

    # Failed Code
    # # Creates the combination list for each letter 
    # CombinationList = []
    # for amountOfCombination in aminoAcid.upper():
    #     InitialList = []
    #     for letter in protein_letters:
    #         if amountOfCombination == letter:
    #             currentCount = 0
    #             actualCount = 0
    #             while actualCount < TotalCombination:
    #                 InitialList.append(protein_letters[letter][currentCount])
    #                 currentCount = currentCount + 1
    #                 if currentCount >= len(protein_letters[letter]):
    #                     currentCount = 0
    #                 actualCount = actualCount + 1 
    #     CombinationList.append(InitialList)

    # # Combines the combination list for each letter to create all the possible combinations
    # CodonList = []
    # for totalCodon in range(0, len(CombinationList[0])):
    #     codonText = []
    #     for codonListCount in range(0, len(CombinationList)):
    #         codonText.append(CombinationList[codonListCount][totalCodon])
    #     CodonList.append(codonText)
    
    # Working Code
    # Takes all the list from the table that needs to be inputted
    CodonList = []
    for amountOfCombination in aminoAcid.upper():
            for letter in protein_letters:
                if amountOfCombination == letter:
                    CodonList.append(protein_letters[amountOfCombination])

    # 
    while len(CodonList) > 1:
        NewCodonList = []
        for codonRow in CodonList[0]:
            for codonRow2 in CodonList[1]:
                SemiText = codonRow + codonRow2
                NewCodonList.append(SemiText)
        CodonList.pop(0)
        CodonList.pop(0)
        CodonList.insert(0, NewCodonList)

    RealCodonList = []
    fullAcid = ""
    acidCount = 0
    for AcidList in CodonList[0]:
        SemiList = []
        for currentAcid in AcidList:
            fullAcid = fullAcid + currentAcid
            acidCount = acidCount + 1
            if acidCount == 3:
                SemiList.append(fullAcid)
                acidCount = 0
                fullAcid = ""
        RealCodonList.append(SemiList)

    # Prints the final result 
    print("Input Aminoacid = " + aminoAcid)
    for codon in RealCodonList:
        codonList = []
        codonFreq = []
        codonName = ""
        for codonTotal in codon:
            codonName = codonName + codonTotal
            if codonTotal not in codonList:
                codonList.append(codonTotal)
                codonFreq.append(1)
            else:
                codonPosition = 0
                for x in codonList:
                    if codonTotal == x:
                        codonFreq[codonPosition] = codonFreq[codonPosition] + 1
                    codonPosition = codonPosition + 1
        print("")
        print("mRNA = " + codonName)
        for codonShow in range(0, len(codonFreq)):
            print(codonList[codonShow] + " = " + str(codonFreq[codonShow]))

# User Interface
def UserInterfaceDNA():
    while True:
        print("")
        print("Please Select which Function to Use: ")
        print("Translate DNA (1)")
        print("Find Codon Frequency (2)")
        print("Exit (Q / 3)")
        inputValue = input("Selection : ")

        if inputValue == "1":
            print("")
            print("Rules: ")
            print("1. String length must be divisible by 3")
            print("2. Most only contain valid dna codon")
            dna_seq = input("Enter DNA Sequence : ")
            Translate_DNA(dna_seq)
        elif inputValue == "2":
            print("")
            print("Rules: ")
            print("1. Contain only valid amino acids")
            aminoAcid = input("Enter Amino Acid : ")
            frequencyCodon(aminoAcid)
        elif inputValue == "Q" or inputValue == "3" or inputValue == "Quit":
            print("Thank you :)")
            return 0
        else:
            print("Please Try Again, Input Value cannot be understood")

# Starts the program
UserInterfaceDNA()