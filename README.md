Identification and determination of the parent of origin for de novo mutations

___________________________________________________________________________________________________________________________________________________
  ![denovo_info_graph](https://github.com/user-attachments/assets/467fad2e-b92b-4f6f-b9c2-ea47e91a2027)


Aims :
de novo mutations are genetic alterations in an organism genome that are not inherited from parents. These mutations are prominent vital for genetic diversity but also lead to disorders.                                               
                          																											            

❖ Vital to profile de novo mutations:

  • Clinical contexts   

  • Understanding population dynamics

  • Understand basic biological processes such as DNA repair

❖ Delineating the parent-of-origin(PoO) of de novo mutations reveals parental specific effects on germline mutation rates.

Inclusion of trio datasets in large scale population studies now allows identifying de novo mutations in the whole genome context. We present, a tool to identify and phase de novo variants from WGS datasets for STRs also SNVs & Indels.

Tool Workflow:

![workflow](https://github.com/user-attachments/assets/a05e9ec0-80f8-4467-951b-0c04b48623ba)


Mutation Model:

![lm](https://github.com/user-attachments/assets/ab43c24b-617d-43bc-a9b0-3b6722b0aec5)



Tool starts with an input joined vcf files of the trio samples. For STRs we currently support outputs of both GangSTR and HipSTR tools. The first step of de novo variant identification picks an allele in the offspring which is not possible to be inherited from parents. We filter out allele drop out cases by checking for reads with the  de novo allele in parent samples. Correspondingly calculate de novo likelihoods using genotype probabilities in parental vs offspring samples. Phasing of  de novo variants is based on two methods, allele sharing and read tracing. In allele sharing we infer the PoO of the de novo allele where the inherited allele has only a single possibility of origin. In read tracing, PoO is determined by building a haplotype with proximal phased heterozygous SNV using read/read-pairs. 




