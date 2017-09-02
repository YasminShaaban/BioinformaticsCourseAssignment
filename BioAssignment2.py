

#-------transfer dna to rna--------
def DNA_to_RNA(dna):
	return dna.replace('T', 'U')

#-------transfer rna to dna--------
def RNA_to_DNA(rna):
	return rna.replace('U', 'T')

#---------------get reverse complement----------------
complementstring= { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',' ':' '}
def reverse_complement(kmer):
	r = ''
	for base in kmer:
		r = complementstring[base] + r
	return r

#---------------complement function------------
def complement(dna):
    basecomplement = {'A': 'T', 'C': 'G','G':'C','g':'c','T':'A','t':'a','N':'N', 'a':'t', 'c':'g',  'n':'n'}
    letters=list(dna)
    letters=[basecomplement[base] for base in basecomplement]
    return ''.join(letters)
#------------------protein translation dictionary-----------

proteintranslation={'A':['GCA','GCC','GCG','GCU'],'S':['UCA','UCC','UCG','UCU'],'L':['UUA','UUG','CUA','CUC','CUG','CUU'],'G':['GGA','GGC','GGG',
                    'GGU'],'P':['CCA','CCG','CCU','CCC'],'R':['CGA','CGC','CGG','CGU','AGA','AGG'],'V':['GUA','GUC','GUG','GUU'],'Stop':['UGA' ,'UAA ','UAG'],
                     'C':['UGC','UGU'],'Q':['CAA','CAG'],'W':['UGG'],'D':['GAC','GAU'],'K':['AAG','AAA'],'T':['ACA','ACC','ACG','ACU'],'I':['AUC','AUA','AUU'],'M':['AUG'],
                     'N':['AAC','AAU'],'F':['UUC','UUU'],'H':['CAC','CAG','CAU'],'E':['GAA'],'Y':['UAU','UAC']}
#'UGA':'Stop' ,'UAA ':'Stop','UAG':'Stop'


#----------------get input-----------------------------------
dnasequence=raw_input("Enter DNA sequence\n")
#print(dnasequence)
proteinsequence=raw_input("Enter Protein Sequence:\n")
#print(proteinsequence)


# --------------------------------------reading 3 frames--------------------------
for n in range (0,3):
    # ---------------------------------ordinary part------------------------------------------
    output=''
#translating to rna
    rnasequence=DNA_to_RNA(dnasequence)

    partpattern=[]

    for i in range (0,len(proteinsequence)):
        partpattern.append('')

    counter=[]
    for i in range (0,len(proteinsequence)):
       counter.append(0)


    for i in range (n,len(rnasequence),1):
      pattern=rnasequence[i:i+3*len(proteinsequence)]
      totalcount=0
      p=0
      for j in range (0,len(pattern),3):
        if (p<len(proteinsequence)):
          counter[p]= 0
          partpattern[p]=pattern[j:j+3]
          amino=proteinsequence[p]
          listamino=proteintranslation[amino]
          for k in range (0,len(listamino),1):
                 if (partpattern[p]==listamino[k]):
                  counter[p]=counter[p]+1
          p+=1

      for ci in range (0,len(proteinsequence)):
          if (counter[ci]==1):
              totalcount+=1
      if (totalcount==len(proteinsequence)):
       output =output + pattern + "  "+"index: "+str(i) +"\n"

    print("frame " + str(n + 1) + ":")
    print(RNA_to_DNA(output))
    #--------------------reverse complement--------------------

    dna = RNA_to_DNA(rnasequence)
    reversecomplement = reverse_complement(dna)

    #print ("reverse complement is" + reversecomplement)

    reversecomplement = DNA_to_RNA(reversecomplement)
    # print(reversecomplement)
    outputrc = ''
    for i in range (0,len(reversecomplement),1):
      pattern=reversecomplement[i:i+3*len(proteinsequence)]
      totalcount=0
      p=0
      for j in range (0,len(pattern),3):
        if (p<len(proteinsequence)):
          counter[p]= 0
          partpattern[p]=pattern[j:j+3]
          amino=proteinsequence[p]
          listamino=proteintranslation[amino]
          for k in range (0,len(listamino),1):
                 if (partpattern[p]==listamino[k]):
                  counter[p]=counter[p]+1
          p+=1

      for ci in range (0,len(proteinsequence)):
          if (counter[ci]==1):
              totalcount+=1
      if (totalcount==len(proteinsequence)):
          for index in range(0, len(dnasequence), 1):
              pattern2 = dnasequence[index:index + 3 * len(proteinsequence)]
              if (pattern2==reverse_complement(RNA_to_DNA(pattern))):
                  index2=index
                  break
          outputrc =outputrc + reverse_complement(RNA_to_DNA(pattern)) +" "+"index: "+str(index2)+"\n"

    print("frame " + str(n + 1) + ":")
    print(RNA_to_DNA(outputrc))






