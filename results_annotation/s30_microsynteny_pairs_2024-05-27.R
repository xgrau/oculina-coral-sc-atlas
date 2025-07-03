# libraries
library("rtracklayer")
library("stringr")
library("igraph")

#### Definitions ####

ogp_fn = "../data/orthology_Anthozoa_plus/orthogroup_pairs.transcripts.csv" # orthology pairs
gff_fo = "../data/reference/" # where can I find gff files?
sps_lis = c("Ocupat","Ocuarb","Spis","Amil","Nvec","Xesp") # species list
out_fn = "results_microsynteny_plus"
dir.create(out_fn, showWarnings = FALSE)

# parsing options
max_dis = 3 # how far away should we go when looking for pairs?

# functions
find_first_introns = function(gff, exon_mark="CDS") {

    # function to identify first intron coordiantes of each transcript in a gtf    
    # find exons & introns
    gfe = gff[gff$type == exon_mark,]
    gfe = gfe[order(gfe$seqid, gfe$start),]
    # intron after each exon (strand-aware) -> first intron will be here, obviously
    gfe$intron_start = ifelse(
        gfe$strand == "+",
        ifelse(gfe$transcript_id == data.table::shift(gfe$transcript_id,1,type="lead"),gfe$end + 1,NA),
        ifelse(gfe$transcript_id == data.table::shift(gfe$transcript_id,1),data.table::shift(gfe$end,1) + 1,NA))
    gfe$intron_end = ifelse(
        gfe$strand == "+",
        ifelse(gfe$transcript_id == data.table::shift(gfe$transcript_id,1,type="lead"),data.table::shift(gfe$start,1,type="lead") - 1,NA),
        ifelse(gfe$transcript_id == data.table::shift(gfe$transcript_id,1),gfe$start - 1,NA))
    # upstream intron coordinates (strand-aware)
    gfe$previntron_start = ifelse(
        gfe$strand == "+",
        ifelse(gfe$transcript_id == data.table::shift(gfe$transcript_id,1),data.table::shift(gfe$end,1) + 1,NA),
        ifelse(gfe$transcript_id == data.table::shift(gfe$transcript_id,1,type="lead"),gfe$end + 1,NA))
    gfe$previntron_end = ifelse(
        gfe$strand == "+",
        ifelse(gfe$transcript_id == data.table::shift(gfe$transcript_id,1),gfe$start - 1,NA),
        ifelse(gfe$transcript_id == data.table::shift(gfe$transcript_id,1,type="lead"),data.table::shift(gfe$start,1,type="lead") - 1,NA))
    
    
    gfe$has_introns = all_duplicated(gfe$transcript_id) 
    gfe$is_first_intron = gfe$has_introns & is.na(gfe$previntron_start)
    gfe_first_introns = gfe[gfe$is_first_intron,c("transcript_id","intron_start","intron_end")]
    colnames(gfe_first_introns) = c("transcript_id","first_intron_start","first_intron_end")
    
    return(gfe_first_introns)
    
}


retrieve_pair_info = function(gff, dis) {
    
    # find neighbouring gene at distance dis
    gff$transcript_b = ifelse(
        gff$seqid == data.table::shift(gff$seqid,n=-dis),
        data.table::shift(gff$transcript_a,n=-dis),
        NA)
    
    # find its orthology group
    gff$orthogroup_b = ifelse(
        gff$seqid == data.table::shift(gff$seqid,n=-dis),
        data.table::shift(gff$orthogroup_a,n=-dis),
        NA)
    
    # find its coordinates
    gff$start_b = ifelse(
        gff$seqid == data.table::shift(gff$seqid,n=-dis),
        data.table::shift(gff$start_a,n=-dis),
        NA)
    gff$end_b = ifelse(
        gff$seqid == data.table::shift(gff$seqid,n=-dis),
        data.table::shift(gff$end_a,n=-dis),
        NA)
    gff$strand_b = ifelse(
        gff$seqid == data.table::shift(gff$seqid,n=-dis),
        data.table::shift(gff$strand_a,n=-dis),
        NA)
    gff$TSS_b = ifelse(
        gff$seqid == data.table::shift(gff$seqid,n=-dis),
        data.table::shift(gff$TSS_a,n=-dis),
        NA)
    
    # first intron coordinates
    gff$first_intron_start_b = ifelse(
        gff$seqid == data.table::shift(gff$seqid,n=-dis),
        data.table::shift(gff$first_intron_start_a,n=-dis),
        NA)
    gff$first_intron_end_b = ifelse(
        gff$seqid == data.table::shift(gff$seqid,n=-dis),
        data.table::shift(gff$first_intron_start_b,n=-dis),
        NA)
    
    # pair info: arrangement
    gff$is_head_to_head = ifelse(
        gff$strand_a == "-",
        ifelse(gff$strand_b == "+",TRUE,FALSE), 
        FALSE)
    gff$is_tail_to_tail = ifelse(
        gff$strand_a == "+",
        ifelse(gff$strand_b == "-",TRUE,FALSE), 
        FALSE)
    gff$is_collinear = !(gff$is_head_to_head | gff$is_tail_to_tail)
    
    # pair distance: TSS
    gff$TSS_distance = abs(gff$TSS_a - gff$TSS_b)
    
    # pair distance: num genes
    gff$gene_distance = dis
    
    # return   
    return(gff)
    
}


# function to identify elements of a vector that are duplicates (a b b c = F T T F)
all_duplicated = function(vec){
    front <- duplicated(vec)
    back <- duplicated(vec, fromLast = TRUE)
    all_dup <- front + back > 0
    return(all_dup)
}

#### Load orthology ####
# orthology
ogp = read.table(ogp_fn, header = FALSE, sep="\t", col.names = c("gene1","gene2"), stringsAsFactors = FALSE)
# add species info
ogp$species1 = stringr::str_split(ogp$gene1, pattern = "_", simplify = TRUE)[,1]
ogp$species2 = stringr::str_split(ogp$gene2, pattern = "_", simplify = TRUE)[,1]
# ogp$species1 [ ogp$species1 == "Ocupat" ] = "Ocubis"
# ogp$species2 [ ogp$species2 == "Ocupat" ] = "Ocubis"


#### Loop through i and j species ####
for (si in sps_lis) {
    
    # load species i
    gi = rtracklayer::readGFF(sprintf("%s/%s_long.annot.gtf", gff_fo,si))
    gig = gi[gi$type == "transcript",]
    # ensure that GFF is sorted by chromosome and starting position
    gig = gig[order(gig$seqid, gig$start),]
    # find TSS
    gig$TSS = ifelse(
        gig$strand == "+",
        gig$start,
        gig$end)
    
    # find first introns
    gig_first_introns = find_first_introns(gff = gi, exon_mark = "CDS")
    gig = merge(gig, gig_first_introns, by.x="transcript_id", by.y="transcript_id", all.x=TRUE, all.y=FALSE)
    
    for (sj in sps_lis) {
        
        if (si != sj) {
            
            message(sprintf("# Pairs %s-%s", si, sj))
            
            # load species j
            gj = rtracklayer::readGFF(sprintf("%s/%s_long.annot.gtf", gff_fo,sj))
            gjg = gj[gj$type == "transcript",]
            # ensure that GFF is sorted by chromosome and starting position
            gjg = gjg[order(gjg$seqid, gjg$start),]
            # find TSS
            gjg$TSS = ifelse(
                gjg$strand == "+",
                gjg$start,
                gjg$end)
            
            # find first introns
            gjg_first_introns = find_first_introns(gff = gj, exon_mark = "CDS")
            gjg = merge(gjg, gjg_first_introns, by.x="transcript_id", by.y="transcript_id", all.x=TRUE, all.y=FALSE)
            
            # subset orthology between i and j
            # subset ortholog pairs
            ogf = ogp[ogp$species1 %in% c(si,sj) & ogp$species2 %in% c(si,sj) & ogp$species1 != ogp$species2 , c("gene1","gene2")]
            # find connected components
            ogf_n = igraph::graph_from_data_frame(ogf, directed = FALSE)
            ogf_n_components = igraph::components(ogf_n)
            ogs = data.frame(
                orthogroup = ogf_n_components$membership,
                gene = names(ogf_n_components$membership))
            
            # add orthogroups to species i and j annotations
            gio = base::merge(gig,ogs, by.x = "transcript_id", by.y="gene", all.x=TRUE, all.y=FALSE)
            gjo = base::merge(gjg,ogs, by.x = "transcript_id", by.y="gene", all.x=TRUE, all.y=FALSE)
            gio = gio[order(gio$seqid, gio$start),]
            gjo = gjo[order(gjo$seqid, gjo$start),]
            
            # drop useless columns
            gio = base::subset(gio, select=c("seqid","start","end","strand","transcript_id","TSS","orthogroup","first_intron_start","first_intron_end"))
            colnames(gio) = c("seqid","start_a","end_a","strand_a","transcript_a","TSS_a","orthogroup_a","first_intron_start_a","first_intron_end_a")
            gjo = base::subset(gjo, select=c("seqid","start","end","strand","transcript_id","TSS","orthogroup","first_intron_start","first_intron_end"))
            colnames(gjo) = c("seqid","start_a","end_a","strand_a","transcript_a","TSS_a","orthogroup_a","first_intron_start_a","first_intron_end_a")
            
            # examine pairs of orthologs at various distances    
            for (dis in 1:max_dis) {
                
                message(sprintf("# Pairs %s-%s | dis = %i | find...", si, sj, dis))
                if (dis == 1) {
                    # init species i and j at dis=1
                    gip = retrieve_pair_info(gff = gio, dis = dis)
                    gjp = retrieve_pair_info(gff = gjo, dis = dis)
                } else {
                    # for subsequent distances...
                    gip = rbind(gip, retrieve_pair_info(gff = gio, dis = dis))
                    gjp = rbind(gjp, retrieve_pair_info(gff = gjo, dis = dis))
                }
            }
            
            # reorder pairs
            gip = gip[order(gip$seqid, gip$start_a, gip$gene_distance),]
            gjp = gjp[order(gjp$seqid, gjp$start_a, gjp$gene_distance),]
            
            # drop pairs where one gene doesn't have orthologs (won't be aligned anyway)
            gip = gip[!is.na(gip$orthogroup_a) & !is.na(gip$orthogroup_b),]
            gjp = gjp[!is.na(gjp$orthogroup_a) & !is.na(gjp$orthogroup_b),]
            
            # drop pairs where both genes belong to the same orthogroup
            gip = gip[gip$orthogroup_a != gip$orthogroup_b,]
            gjp = gjp[gjp$orthogroup_a != gjp$orthogroup_b,]
            
            # OG pair id (internally sorted, for uniqueness)
            gip$og_pair_id = ifelse(
                gip$orthogroup_a <= gip$orthogroup_b, 
                paste(gip$orthogroup_a , gip$orthogroup_b), 
                paste(gip$orthogroup_b , gip$orthogroup_a))
            gjp$og_pair_id = ifelse(
                gjp$orthogroup_a <= gjp$orthogroup_b, 
                paste(gjp$orthogroup_a , gjp$orthogroup_b), 
                paste(gjp$orthogroup_b , gjp$orthogroup_a))
            
            # merge pairs
            message(sprintf("# Pairs %s-%s | dis up to %i | report...", si, sj, max_dis))
            dpa = base::merge(gip, gjp, by.x="og_pair_id", by.y="og_pair_id", all.x=FALSE, all.y=FALSE, suffixes=c("_i","_j"))
            
            # double pair info
            # is this a unique pair?
            dpa$is_unique_pair = !all_duplicated(dpa$og_pair_id)
            
            # report
            write.table(dpa, file = sprintf("%s/pairs_%s-%s.csv", out_fn, si, sj), quote = FALSE, row.names = FALSE, sep="\t")
            
        }
    }
}