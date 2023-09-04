#Created by Siwen He (siwenhe@cqu.edu.cn) at 2023/09/04
#We used the following R scripts that were adapted from previous work relating beta diversity partitioning (e.g., see a case study in Khattar et al. 2021, Ecography, 44, 1391-1402) to partition total beta diversity into its S scale components

#Function developed to Parition BDtotal into S scale components

#We measured S as: 
##S = 2n+1, c = 1; if presence of a spatial category for the sampling sites (Category=1) 
## OR S = n+1, c = 0; if absence of a spatial category for the sampling sites (Category=0) 
### n are the numbers of nested hierarchies (Nested_hierarchy) included in the study design 
### For simplicity, we included only the first four nested hierarchies. Additionally, we believe that most studies would not include more than four nested hierarchies in their study design  
### However,one can expand it to other nested hierarchies (Nested_hierarchy = 5,6,...) by modifing our R scripts 

#Please note that we used stream network comprised sites (communities) with a spatial category (i.e., Strahler stream order, c = 1) or lack a spatial category (c = 0) and ecoregions with one (i.e. LevelIV, n =1) to four nested hierarchies (i.e., Level I, Level II, Level III and Level IV, Omernik and Griffith 2014, Environ. Manage. 54:1249-1266, n = 4) as specific examples 
##So it requires from the user a very specific type of naming for communities
###Specified name the communities (row names) for instance as follows: LevelI_LevelII_LevelIII_LevelIV_Strahler stream order,if category = 1 and Nested_hierarchy = 4

#### Category=1 (presence of a spatial category)
#### Nested_hierarchy = 1 (Using LevelIV as an example, Parition BDtotal into BDam.LIV,BDWA.LIV and BDwi.LIV)
#### Nested_hierarchy = 2 (Using LevelIII and LevelIV as an example,Parition BDtotal into BDam.LIII, BDam.LIV, BDWA.LIII, BDWA.LIV and BDwi.LIV)
#### Nested_hierarchy = 3 (Using LevelII, LevelIII and LevelIV as an examplet,Parition BDtotal into BDam.LII, BDam.LIII, BDam.LIV, BDWA.LII, BDWA.LIII, BDWA.LIV and BDwi.LIV)
#### Nested_hierarchy = 4 (Using LevelI, LevelII, LevelIII and LevelIV as an example, Parition BDtotal into BDam.LI, BDam.LII, BDam.LIII, BDam.LIV, BDWA.LI, BDWA.LII, BDWA.LIII, BDWA.LIV and BDwi.LIV)

#### Category=0 (absence of a spatial category)
#### Nested_hierarchy = 1 (Using LevelIV as an example, Parition BDtotal into BDAm.LIV and BDwi.LIV)
#### Nested_hierarchy = 2 (Using LevelIII and LevelIV as an example,Parition BDtotal into BDAm.LIII, BDAm.LIV and BDwi.LIV)
#### Nested_hierarchy = 3 (Using LevelII, LevelIII and LevelIV as an examplet,Parition BDtotal into BDAm.LII, BDAm.LIII, BDAm.LIV and BDwi.LIV)
#### Nested_hierarchy = 4 (Using LevelI, LevelII, LevelIII and LevelIV as an example, Parition BDtotal into BDAm.LI, BDAm.LII, BDAm.LIII, BDAm.LIV and BDwi.LIV)
#### Please note that the “BDAm.LI, BDAm.LII, BDAm.LIII, BDAm.LIV”(indicates purely “among” plus “within and among”) used in the Category=0 is different to the “BDam.LI, BDam.LII, BDam.LIII, BDam.LIV” (indicates purely “among”) used in the Category=1

# It retunrs a list where:
# BD_Partition = Partition of BDtotal across spatial scales 

#Total_entries_per_scale = Number of entries in D representing beta diversity in each hierarchical scale

#x= Species-by-Strahler-order-by-ecoregion matrix, if category=1
#OR x= Species-by-Site-by-ecoregion, if category=0
#Sample names (row names) indicate the hierarchical ecoregional (LI, LII, LIII, LIV) and network (Strahler-order, SO) location of samples as follows: LI_LII_LIII_LIV_SO,if category=1
#OR Sample names (row names) indicate the hierarchical ecoregional (LI, LII, LIII, LIV) and site as follows: LI_LII_LIII_LIV_SI,if category=0



library(adespatial) #Load required library #requires function "beta.div.comp"
#Component = "Total" Calculates Sorensen Dissimilarity (Baselga family); "Turnover" Calculates Turnover; "Nestedness" Caculates Nestdedness
#Category=1 presence of a spatial category; Category=0 absence of a spatial category
#Nested_hierarchy = 1 one nested hierarchies ; Nested_hierarchy = 2 two nested hierarchies ;Nested_hierarchy = 3 three nested hierarchies ; Nested_hierarchy = 4 four nested hierarchies


BD_partition<-function(x,Component,Category,Nested_hierarchy){
  
  if(Component=="Total"){
    Diss_matrix<-beta.div.comp(x,coef="BS",quant=T)$D
    Diss_total<-data.frame(t(combn(as.character(rownames(x)),2)), as.numeric(Diss_matrix))
    names(Diss_total) <- c("c1", "c2", "distance")
  }
  
  if(Component=="Turnover"){
    Diss_matrix<-beta.div.comp(x,coef="BS",quant=T)$repl
    Diss_total<-data.frame(t(combn(as.character(rownames(x)),2)), as.numeric(Diss_matrix))
    names(Diss_total) <- c("c1", "c2", "distance")
  }
  
  if(Component=="Nestedness"){
    Diss_matrix<-beta.div.comp(x,coef="BS",quant=T)$rich
    Diss_total<-data.frame(t(combn(as.character(rownames(x)),2)), as.numeric(Diss_matrix))
    names(Diss_total) <- c("c1", "c2", "distance")
  }
  
  # presence of a spatial category
  if(Category==1){
    
    # Naming the scale of each entry in Diss_matrix
    # one nested_hierarchy (LevelIV) included in the study design
    if(Nested_hierarchy==1){
      Dim<-vector()
      for(j in 1:nrow(Diss_total)){
        s1.letter = strsplit(as.character(Diss_total[j,1]), split = "_")[[1]]
        s2.letter = strsplit(as.character(Diss_total[j,2]), split= "_") [[1]]
        if(s1.letter[1] == s2.letter[1] & s1.letter[2] != s2.letter[2]){
          Dim[i]<- "wi.LIV"
        }else{
          if(s1.letter[1] != s2.letter[1] & s1.letter[2] == s2.letter[2]){
            Dim[i]<- "am.LIV"
          }else{
            Dim[i]<- "WA.LIV"
          }}
      }
      Diss_total$Dim<-Dim
      
      Total_entries_per_scale<-table(Dim)# Total entries in D representing each scale
      
      BDwi.LIV<-sum(subset(Diss_total,subset=Dim=="wi.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LIV<-sum(subset(Diss_total,subset=Dim=="am.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDWA.LIV<-sum(subset(Diss_total,subset=Dim=="WA.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      
      BD_per_scale<- c(BDam.LIV, BDWA.LIV,BDwi.LIV)
      names(BD_per_scale)<-c("am.LIV",  "WA.LIV",  "wi.LIV" )
      
      BD_partition<-BD_per_scale/Total_entries_per_scale#Avg. Cont to BDtotal(i.e., BD_per_scale/Total entries per dimension)
      names(BD_partition)<-c("BDam.LIV","BDWA.LIV", "BDwi.LIV")
      
      Result<-list(BD_partition,Total_entries_per_scale)
      names(Result)<-c("BD_partition","Total entries in D per scale")
      return(Result) 
    }
    
    # Naming the scale of each entry in Diss_matrix
    # two nested_hierarchies (LevelIII and LevelIV) included in the study design
    if(Nested_hierarchy==2){
      Dim<-vector()
      for(j in 1:nrow(Diss_total)){
        s1.letter = strsplit(as.character(Diss_total[j,1]), split = "_")[[1]]
        s2.letter = strsplit(as.character(Diss_total[j,2]), split= "_") [[1]]
        if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]){
          Dim[j]<- "wi.LIV"
        }else{
          if(s1.letter[1] == s2.letter[1] & s1.letter[2] != s2.letter[2]& s1.letter[3] == s2.letter[3]){
            Dim[j]<- "am.LIV"
          }else{ 
            if(s1.letter[1] == s2.letter[1] & s1.letter[2] != s2.letter[2]& s1.letter[3] != s2.letter[3]){
              Dim[j]<- "WA.LIV"
            }else{
              if(s1.letter[1] != s2.letter[1]& s1.letter[3] == s2.letter[3]){
                Dim[j]<- "am.LIII"
              }else{
                 Dim[j]<- "WA.LIII"
                      }}}}
      }
      Diss_total$Dim<-Dim
      
      Total_entries_per_scale<-table(Dim)# Total entries in D representing each scale
      
      BDwi.LIV<-sum(subset(Diss_total,subset=Dim=="wi.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LIV<-sum(subset(Diss_total,subset=Dim=="am.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDWA.LIV<-sum(subset(Diss_total,subset=Dim=="WA.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LIII<-sum(subset(Diss_total,subset=Dim=="am.LIII")[,3])/(nrow(x)*(nrow(x)-1))
      BDWA.LIII<-sum(subset(Diss_total,subset=Dim=="WA.LIII")[,3])/(nrow(x)*(nrow(x)-1))
     
      
      BD_per_scale<- c(BDam.LIII,BDam.LIV,BDWA.LIII, BDWA.LIV,BDwi.LIV)
      names(BD_per_scale)<-c("am.LIII","am.LIV","WA.LIII",  "WA.LIV",  "wi.LIV" )
      
      BD_partition<-BD_per_scale/Total_entries_per_scale#Avg. Cont to BDtotal(i.e., BD_per_scale/Total entries per dimension)
      names(BD_partition)<-c("BDam.LIII","BDam.LIV","BDWA.LIII", "BDWA.LIV", "BDwi.LIV")
      
      Result<-list(BD_partition,Total_entries_per_scale)
      names(Result)<-c("BD_partition","Total entries in D per scale")
      return(Result) 
    }
    
    # Naming the scale of each entry in Diss_matrix
    # three nested_hierarchies (LevelII, LevelIII and LevelIV) included in the study design
    if(Nested_hierarchy==3){
      Dim<-vector()
      for(j in 1:nrow(Diss_total)){
        s1.letter = strsplit(as.character(Diss_total[j,1]), split = "_")[[1]]
        s2.letter = strsplit(as.character(Diss_total[j,2]), split= "_") [[1]]
        if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] == s2.letter[3]){
          Dim[j]<- "wi.LIV"
        }else{
          if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] != s2.letter[3]& s1.letter[4] == s2.letter[4]){
            Dim[j]<- "am.LIV"
          }else{ 
            if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] != s2.letter[3]& s1.letter[4] != s2.letter[4]){
              Dim[j]<- "WA.LIV"
            }else{
              if(s1.letter[1] == s2.letter[1] & s1.letter[2] != s2.letter[2]& s1.letter[4] == s2.letter[4]){
                Dim[j]<- "am.LIII"
              }else{
                if(s1.letter[1] == s2.letter[1] & s1.letter[2] != s2.letter[2]& s1.letter[4] != s2.letter[4]){
                  Dim[j]<- "WA.LIII"
                }else{
                  if(s1.letter[1] != s2.letter[1] & s1.letter[4] == s2.letter[4]){
                    Dim[j]<- "am.LII"
                  }else{
                    Dim[j]<- "WA.LII"
                      }}}}}}
      }
      Diss_total$Dim<-Dim
      
      Total_entries_per_scale<-table(Dim)# Total entries in D representing each scale
      
      BDwi.LIV<-sum(subset(Diss_total,subset=Dim=="wi.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LIV<-sum(subset(Diss_total,subset=Dim=="am.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDWA.LIV<-sum(subset(Diss_total,subset=Dim=="WA.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LIII<-sum(subset(Diss_total,subset=Dim=="am.LIII")[,3])/(nrow(x)*(nrow(x)-1))
      BDWA.LIII<-sum(subset(Diss_total,subset=Dim=="WA.LIII")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LII<-sum(subset(Diss_total,subset=Dim=="am.LII")[,3])/(nrow(x)*(nrow(x)-1))
      BDWA.LII<-sum(subset(Diss_total,subset=Dim=="WA.LII")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LI<-sum(subset(Diss_total,subset=Dim=="am.LI")[,3])/(nrow(x)*(nrow(x)-1))
      BDWA.LI<-sum(subset(Diss_total,subset=Dim=="WA.LI")[,3])/(nrow(x)*(nrow(x)-1))
      
      BD_per_scale<- c(BDam.LII,BDam.LIII,BDam.LIV,BDWA.LII,BDWA.LIII, BDWA.LIV,BDwi.LIV)
      names(BD_per_scale)<-c("am.LII","am.LIII","am.LIV","WA.LII","WA.LIII",  "WA.LIV",  "wi.LIV" )
      
      BD_partition<-BD_per_scale/Total_entries_per_scale#Avg. Cont to BDtotal(i.e., BD_per_scale/Total entries per dimension)
      names(BD_partition)<-c("BDam.LII","BDam.LII","BDam.LIV",
                             "BDWA.LII","BDWA.LIII", "BDWA.LIV", "BDwi.LIV")
      
      Result<-list(BD_partition,Total_entries_per_scale)
      names(Result)<-c("BD_partition","Total entries in D per scale")
      return(Result) 
    }

    # Naming the scale of each entry in Diss_matrix
    # four nested_hierarchies (LevelI, LevelII, LevelIII and LevelIV) included in the study design
    if(Nested_hierarchy==4){
        Dim<-vector()
        for(j in 1:nrow(Diss_total)){
          s1.letter = strsplit(as.character(Diss_total[j,1]), split = "_")[[1]]
          s2.letter = strsplit(as.character(Diss_total[j,2]), split= "_") [[1]]
          if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] == s2.letter[3]& s1.letter[4] == s2.letter[4]){
            Dim[j]<- "wi.LIV"
          }else{
            if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] == s2.letter[3]& s1.letter[4] != s2.letter[4]& s1.letter[5] == s2.letter[5]){
              Dim[j]<- "am.LIV"
            }else{ 
              if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] == s2.letter[3]& s1.letter[4] != s2.letter[4]& s1.letter[5] != s2.letter[5]){
                Dim[j]<- "WA.LIV"
              }else{
                if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] != s2.letter[3]& s1.letter[5] == s2.letter[5]){
                  Dim[j]<- "am.LIII"
                }else{
                  if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] != s2.letter[3]& s1.letter[5] != s2.letter[5]){
                    Dim[j]<- "WA.LIII"
                  }else{
                    if(s1.letter[1] == s2.letter[1] & s1.letter[2] != s2.letter[2]& s1.letter[5] == s2.letter[5]){
                      Dim[j]<- "am.LII"
                    }else{
                      if(s1.letter[1] == s2.letter[1] & s1.letter[2] != s2.letter[2]& s1.letter[5] != s2.letter[5]){
                        Dim[j]<- "WA.LII"
                      }else{
                        if(s1.letter[1] != s2.letter[1] & s1.letter[5] == s2.letter[5]){
                          Dim[j]<- "am.LI"
                        }else{
                          Dim[j]<- "WA.LI"
                        }}}}}}}}
        }
        Diss_total$Dim<-Dim
        
        Total_entries_per_scale<-table(Dim)# Total entries in D representing each scale
        
        BDwi.LIV<-sum(subset(Diss_total,subset=Dim=="wi.LIV")[,3])/(nrow(x)*(nrow(x)-1))
        BDam.LIV<-sum(subset(Diss_total,subset=Dim=="am.LIV")[,3])/(nrow(x)*(nrow(x)-1))
        BDWA.LIV<-sum(subset(Diss_total,subset=Dim=="WA.LIV")[,3])/(nrow(x)*(nrow(x)-1))
        BDam.LIII<-sum(subset(Diss_total,subset=Dim=="am.LIII")[,3])/(nrow(x)*(nrow(x)-1))
        BDWA.LIII<-sum(subset(Diss_total,subset=Dim=="WA.LIII")[,3])/(nrow(x)*(nrow(x)-1))
        BDam.LII<-sum(subset(Diss_total,subset=Dim=="am.LII")[,3])/(nrow(x)*(nrow(x)-1))
        BDWA.LII<-sum(subset(Diss_total,subset=Dim=="WA.LII")[,3])/(nrow(x)*(nrow(x)-1))
        BDam.LI<-sum(subset(Diss_total,subset=Dim=="am.LI")[,3])/(nrow(x)*(nrow(x)-1))
        BDWA.LI<-sum(subset(Diss_total,subset=Dim=="WA.LI")[,3])/(nrow(x)*(nrow(x)-1))
        
        BD_per_scale<- c(BDam.LI, BDam.LII,BDam.LIII,BDam.LIV,BDWA.LI,BDWA.LII,BDWA.LIII, BDWA.LIV,BDwi.LIV)
        names(BD_per_scale)<-c("am.LI","am.LII","am.LIII","am.LIV","WA.LI","WA.LII","WA.LIII",  "WA.LIV",  "wi.LIV" )
        
        BD_partition<-BD_per_scale/Total_entries_per_scale#Avg. Cont to BDtotal(i.e., BD_per_scale/Total entries per dimension)
        names(BD_partition)<-c("BDam.LI","BDam.LII","BDam.LII","BDam.LIV","BDWA.LI",
                               "BDWA.LII","BDWA.LIII", "BDWA.LIV", "BDwi.LIV")
        
        Result<-list(BD_partition,Total_entries_per_scale)
        names(Result)<-c("BD_partition","Total entries in D per scale")
        return(Result) 
    }
  }
  
  
  # absence of a spatial category
  if(Category==0){
    
    # Naming the scale of each entry in Diss_matrix
    # one nested_hierarchy (LevelIV) included in the study design
    if(Nested_hierarchy==1){
      Dim<-vector()
      for(j in 1:nrow(Diss_total)){
        s1.letter = strsplit(as.character(Diss_total[j,1]), split = "_")[[1]]
        s2.letter = strsplit(as.character(Diss_total[j,2]), split= "_") [[1]]
        if(s1.letter[1] == s2.letter[1]){
          Dim[i]<- "wi.LIV"
        }else{
           Dim[i]<- "Am.LIV"
          }
      }
      Diss_total$Dim<-Dim
      
      Total_entries_per_scale<-table(Dim)# Total entries in D representing each scale
      
      BDwi.LIV<-sum(subset(Diss_total,subset=Dim=="wi.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LIV<-sum(subset(Diss_total,subset=Dim=="Am.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      
      
      BD_per_scale<- c(BDam.LIV, BDwi.LIV)
      names(BD_per_scale)<-c("Am.LIV", "wi.LIV" )
      
      BD_partition<-BD_per_scale/Total_entries_per_scale#Avg. Cont to BDtotal(i.e., BD_per_scale/Total entries per dimension)
      names(BD_partition)<-c("BDAm.LIV", "BDwi.LIV")
      
      Result<-list(BD_partition,Total_entries_per_scale)
      names(Result)<-c("BD_partition","Total entries in D per scale")
      return(Result) 
    }
    
    # Naming the scale of each entry in Diss_matrix
    # two nested_hierarchies (LevelIII and LevelIV) included in the study design
    if(Nested_hierarchy==2){
      Dim<-vector()
      for(j in 1:nrow(Diss_total)){
        s1.letter = strsplit(as.character(Diss_total[j,1]), split = "_")[[1]]
        s2.letter = strsplit(as.character(Diss_total[j,2]), split= "_") [[1]]
        if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]){
          Dim[j]<- "wi.LIV"
        }else{
          if(s1.letter[1] == s2.letter[1] & s1.letter[2] != s2.letter[2]){
            Dim[j]<- "Am.LIV"
          }else{ 
            Dim[j]<- "Am.LIII"
              }}
      }
      Diss_total$Dim<-Dim
      
      Total_entries_per_scale<-table(Dim)# Total entries in D representing each scale
      
      BDwi.LIV<-sum(subset(Diss_total,subset=Dim=="wi.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LIV<-sum(subset(Diss_total,subset=Dim=="Am.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LIII<-sum(subset(Diss_total,subset=Dim=="Am.LIII")[,3])/(nrow(x)*(nrow(x)-1))
      
      
      BD_per_scale<- c(BDam.LIII,BDam.LIV,BDwi.LIV)
      names(BD_per_scale)<-c("Am.LIII","Am.LIV", "wi.LIV" )
      
      BD_partition<-BD_per_scale/Total_entries_per_scale#Avg. Cont to BDtotal(i.e., BD_per_scale/Total entries per dimension)
      names(BD_partition)<-c("BDAm.LIII","BDAm.LIV", "BDwi.LIV")
      
      Result<-list(BD_partition,Total_entries_per_scale)
      names(Result)<-c("BD_partition","Total entries in D per scale")
      return(Result) 
    }
    
    # Naming the scale of each entry in Diss_matrix
    # three nested_hierarchies (LevelII, LevelIII and LevelIV) included in the study design
    if(Nested_hierarchy==3){
      Dim<-vector()
      for(j in 1:nrow(Diss_total)){
        s1.letter = strsplit(as.character(Diss_total[j,1]), split = "_")[[1]]
        s2.letter = strsplit(as.character(Diss_total[j,2]), split= "_") [[1]]
        if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] == s2.letter[3]){
          Dim[j]<- "wi.LIV"
        }else{
          if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] != s2.letter[3]){
            Dim[j]<- "Am.LIV"
          }else{ 
            if(s1.letter[1] == s2.letter[1] & s1.letter[2] != s2.letter[2]){
              Dim[j]<- "Am.LIII"
            }else{
              Dim[j]<- "Am.LII"
                  }}}
      }
      Diss_total$Dim<-Dim
      
      Total_entries_per_scale<-table(Dim)# Total entries in D representing each scale
      
      BDwi.LIV<-sum(subset(Diss_total,subset=Dim=="wi.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LIV<-sum(subset(Diss_total,subset=Dim=="Am.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LIII<-sum(subset(Diss_total,subset=Dim=="Am.LIII")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LII<-sum(subset(Diss_total,subset=Dim=="Am.LII")[,3])/(nrow(x)*(nrow(x)-1))
      
      
      BD_per_scale<- c(BDam.LII,BDam.LIII,BDam.LIV,BDwi.LIV)
      names(BD_per_scale)<-c("Am.LII","Am.LIII","Am.LIV", "wi.LIV" )
      
      BD_partition<-BD_per_scale/Total_entries_per_scale#Avg. Cont to BDtotal(i.e., BD_per_scale/Total entries per dimension)
      names(BD_partition)<-c("BDAm.LII","BDAm.LII","BDAm.LIV", "BDwi.LIV")
      
      Result<-list(BD_partition,Total_entries_per_scale)
      names(Result)<-c("BD_partition","Total entries in D per scale")
      return(Result) 
    }
    
    # Naming the scale of each entry in Diss_matrix
    # four nested_hierarchies (LevelI, LevelII, LevelIII and LevelIV) included in the study design
    if(Nested_hierarchy==4){
      Dim<-vector()
      for(j in 1:nrow(Diss_total)){
        s1.letter = strsplit(as.character(Diss_total[j,1]), split = "_")[[1]]
        s2.letter = strsplit(as.character(Diss_total[j,2]), split= "_") [[1]]
        if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] == s2.letter[3]& s1.letter[4] == s2.letter[4]){
          Dim[j]<- "wi.LIV"
        }else{
          if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] == s2.letter[3]& s1.letter[4] != s2.letter[4]){
            Dim[j]<- "Am.LIV"
          }else{ 
            if(s1.letter[1] == s2.letter[1] & s1.letter[2] == s2.letter[2]& s1.letter[3] != s2.letter[3]){
              Dim[j]<- "Am.LIII"
            }else{
              if(s1.letter[1] == s2.letter[1] & s1.letter[2] != s2.letter[2]){
                Dim[j]<- "Am.LII"
              }else{
                Dim[j]<- "Am.LI"
                      }}}}
      }
      Diss_total$Dim<-Dim
      
      Total_entries_per_scale<-table(Dim)# Total entries in D representing each scale
      
      BDwi.LIV<-sum(subset(Diss_total,subset=Dim=="wi.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LIV<-sum(subset(Diss_total,subset=Dim=="Am.LIV")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LIII<-sum(subset(Diss_total,subset=Dim=="Am.LIII")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LII<-sum(subset(Diss_total,subset=Dim=="Am.LII")[,3])/(nrow(x)*(nrow(x)-1))
      BDam.LI<-sum(subset(Diss_total,subset=Dim=="Am.LI")[,3])/(nrow(x)*(nrow(x)-1))
      
      
      BD_per_scale<- c(BDam.LI, BDam.LII,BDam.LIII,BDam.LIV,BDwi.LIV)
      names(BD_per_scale)<-c("Am.LI","Am.LII","Am.LIII","Am.LIV", "wi.LIV" )
      
      BD_partition<-BD_per_scale/Total_entries_per_scale#Avg. Cont to BDtotal(i.e., BD_per_scale/Total entries per dimension)
      names(BD_partition)<-c("BDAm.LI","BDAm.LII","BDAm.LII","BDAm.LIV","BDwi.LIV")
      
      Result<-list(BD_partition,Total_entries_per_scale)
      names(Result)<-c("BD_partition","Total entries in D per scale")
      return(Result) 
    }
  }
}

#### Examples#####

Samples<-32 
# 2 Strahler orders (SO1, SO2) in each of 16 LIV ecoregions, 
#belong to 8 LIII ecoregions nested in 4 LII ecoregions on 2 LI ecoregions (LI1, LI2)
LI<-c(rep("LI1",times=16),rep("LI2",times=16))
LI_LII<-c(rep("LI1.LII1",times=8),rep("LI1.LII2",times=8),
          rep("LI2.LII1",times=8),rep("LI2.LII2",times=8))
LI_LII_LIII<-c(rep("LI1.LII1.LIII1",times=4),rep("LI1.LII1.LIII2",times=4),
               rep("LI1.LII2.LIII1",times=4),rep("LI1.LII2.LIII2",times=4),
               rep("LI2.LII1.LIII1",times=4),rep("LI2.LII1.LIII2",times=4),
               rep("LI2.LII2.LIII1",times=4),rep("LI2.LII2.LIII2",times=4))
LI_LII_LIII_LIV<-c(rep("LI1.LII1.LIII1.LIV1",times=2),
                   rep("LI1.LII1.LIII1.LIV2",times=2),
                   rep("LI1.LII1.LIII2.LIV1",times=2),
                   rep("LI1.LII1.LIII2.LIV2",times=2),
                   rep("LI1.LII2.LIII1.LIV1",times=2),
                   rep("LI1.LII2.LIII1.LIV2",times=2),
                   rep("LI1.LII2.LIII2.LIV1",times=2),
                   rep("LI1.LII2.LIII2.LIV2",times=2),
                   rep("LI2.LII1.LIII1.LIV1",times=2),
                   rep("LI2.LII1.LIII1.LIV2",times=2),
                   rep("LI2.LII1.LIII2.LIV1",times=2),
                   rep("LI2.LII1.LIII2.LIV2",times=2),
                   rep("LI2.LII2.LIII1.LIV1",times=2),
                   rep("LI2.LII2.LIII1.LIV2",times=2),
                   rep("LI2.LII2.LIII2.LIV1",times=2),
                   rep("LI2.LII2.LIII2.LIV2",times=2))
SO<-c("SO1","SO2","SO1","SO2", "SO1","SO2","SO1","SO2",
      "SO1","SO2","SO1","SO2", "SO1","SO2","SO1","SO2",
      "SO1","SO2","SO1","SO2", "SO1","SO2","SO1","SO2",
      "SO1","SO2","SO1","SO2", "SO1","SO2","SO1","SO2")

library(stringr)#Load required librariy ## requires function "str_c"
# LI_LII_LIII_LIV_SO = Sample names (row names)

LI_LII_LIII_LIV_SO<-c(str_c(LI,"_",LI_LII, "_", LI_LII_LIII,"_", LI_LII_LIII_LIV,
                            "_", SO))

Species<-50
Regional_abundance<-500

x<-matrix(0,nrow = Samples,ncol=Species,
          dimnames = list(LI_LII_LIII_LIV_SO,c(paste0("species",1:Species))))

#Lognormal SAD##Load required librariy 
library(mobsim) #requires function "sim_sad"

lognormal_pool<-sim_sad(s_pool=Species,n_sim=Regional_abundance,sad_type="lnorm",
                        sad_coef = list("meanlog"=5,"sdlog"=0.5))

#Populating x
for(i in 1:length(lognormal_pool)){
  
  Distribution<-(table(sample(rownames(x),as.numeric(lognormal_pool)[i],replace=T)))
  Sites<-names(Distribution)  
  Local_abundances<-as.numeric(Distribution)
  x[Sites,i]<-Local_abundances
  
}

BDtotal_Result<-BD_partition(x,Component = "Total", Category=1,Nested_hierarchy=4)

#end