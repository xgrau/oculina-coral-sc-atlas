# Datasets

## Single-cell atlases

List of species:

```tsv
# Dataset ID	Species name	Data source
Ocupat	Oculina patagonica	This study
Ocuarb	Oculina arbuscula	https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1122932
Spin	Stylophora pistillata	This study; reference data from https://pubmed.ncbi.nlm.nih.gov/33945788/
Amil	Acropora millepora	This study
Nvec	Nematostella vectensis	https://pubmed.ncbi.nlm.nih.gov/29856957/
Xesp	Xenia sp.	https://pubmed.ncbi.nlm.nih.gov/32555454/
```

## Comparative genomics

### Anthozoa species set

List of species:

```tsv
# Species ID	Species name	Taxonomy	Data source
Actieq	Actinia equina	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Actiniaria, Actiniidae	Reef Genomics v1 http://aequ.reefgenomics.org/download/ ; NCBI GCA_011057435.1
Metsen	Metridium senile	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Actiniaria, Acuticulata, Metridiidae	NCBI GCA_949775045.1, jaMetSeni4.1 DTOL; in-home annotation
Dialin	Diadumene lineata	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Actiniaria, Acuticulata, Diadumenidae	NCBI GCA_918843875.1, jaDiaLine6.1 DTOL; in-home annotation
Exapal	Exaiptasia diaphana	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Actiniaria, Acuticulata, Aiptasiidae	NCBI GCF_001417965.1_Aiptasia_genome_1.1
Nvec	Nematostella vectensis	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Actiniaria, Edwardsiidae	chromosomal assembly jaNemVect1.1, GCF_932526225.1; Combined JGI + Vienna annotation (https://doi.org/10.1016/j.cell.2018.05.019)
Scocal	Scolanthus callimorphus	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Actiniaria, Edwardsiidae	Stowers Institute SIMRBASE version Scal100_v1.20200813
Adig	Acropora digitifera	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Acroporidae	NCBI GCF_000222465
Amil	Acropora millepora	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Acroporidae	NCBI assembly: GCA_013753865.1 (v2.1); annotation from Fuller 2020 (via https://drive.google.com/file/d/1ww7rbr6v8676Mul8OxoHkTp3fau0h27F/view?usp=sharing)
Gfas	Galaxea fascicularis	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Euphylliidae	Reef Genomics v1 http://gfas.reefgenomics.org/download/
Fspp	Fungia sp.	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Fungiidae	Reef Genomics v1 http://ffun.reefgenomics.org/download/
Ocupat	Oculina patagonica	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Oculinidae	This study
Ocuarb	Oculina arbuscula	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Oculinidae	NCBI GCA_964656845.1: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964656845.1/
Gasp	Coelastrea aspera	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Merulinidae	Reef Genomics v1 http://gasp.reefgenomics.org/download/; formerly Goniastrea aspera
Spis	Stylophora pistillata	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Pocilloporidae	NCBI GCA_002571385.1 Stylophora pistillata v1 + reefgenomics (see Methods from Levy et al Cell 2021)
Pocdam	Pocillopora damicornis	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Pocilloporidae	NCBI GCF_003704095.1_ASM370409v1
Dgig	Dendronephthya gigantea	Metazoa, Cnidaria, Anthozoa, Octocorallia	NCBI GCF_004324835.1_DenGig_1.0
Xesp	Xenia sp.	Metazoa, Cnidaria, Anthozoa, Octocorallia	Carnegie Coral & Marine Organisms database, https://cmo.carnegiescience.edu/data/
```

Species tree:

```bash
((((Actieq,((Metsen,Dialin)MetrDiad,Exapal)Acuticulata)Enthemonae,(Nvec,Scocal)Edwardsiidae)Actiniaria,(((Adig,Amil)Acroporidae,Gfas)AcroEuphylliidae,(((Fspp,(Ocupat,Ocuarb)Oculinidae)OcuFun,Gasp)FunMerulinidae,(Spis,Pocdam)Pocilloporidae)PocFunMerulinidae)Scleractinia)Hexacorallia,(Dgig,Xesp)Octocorallia)Anthozoa;
```

### Anthozoa species set, species with chromosome-level assemblies

```tsv
# Species ID	Species name	Taxonomy	Data source
Amil	Acropora millepora	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Acroporidae	NCBI assembly: GCA_013753865.1 (v2.1); annotation from Fuller 2020 (via https://drive.google.com/file/d/1ww7rbr6v8676Mul8OxoHkTp3fau0h27F/view?usp=sharing)
Acrcer	Acropora cervicornis	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Acroporidae	NCBI GCA_964034985.1 jaAcrCerv1.1; in-home annotation
Acrpal	Acropora palmata	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Acroporidae	NCBI GCA_964030605.1 jaAcrPala1.1; in-home annotation
Ocupat	Oculina patagonica	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Oculinidae	This study
Ocuarb	Oculina arbuscula	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Oculinidae	NCBI GCA_964656845.1: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964656845.1/
Porlut	Porites lutea	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Poritidae	NCBI GCA_958299795.1, annotated from ERR12708749; in-home annotation
Pocver	Pocillopora verrucosa	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Pocilloporidae	NCBI GCA_036669915.2_ASM3666991v2, annotated from SRR11880672; in-home annotation
Nvec	Nematostella vectensis	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Actiniaria, Edwardsiidae	chromosomal assembly jaNemVect1.1, GCF_932526225.1; Combined JGI + Vienna annotation (https://doi.org/10.1016/j.cell.2018.05.019)
Xesp	Xenia sp.	Metazoa, Cnidaria, Anthozoa, Octocorallia	Carnegie Coral & Marine Organisms database, https://cmo.carnegiescience.edu/data/
Rhoesc	Rhopilema esculentum	Metazoa, Cnidaria, Medusozoa, Scyphozoa, Rhizostomeae	NCBI GCA_013076305.1_ASM1307630v1 / annotation: https://github.com/nongwy/JellyfishGenomeData (scaffolds need to be mapped to NCBI names)
Hvul	Hydra vulgaris	Metazoa, Cnidaria, Medusozoa, Hydrozoa, Anthoathecata, Hydridae	NIH https://research.nhgri.nih.gov/hydra/
```

Species tree:

```bash
((((((Ocupat,Ocuarb)Oculinidae,Pocver),((Amil,(Acrcer,Acrpal)),Porlut)),Nvec),Xesp),(Hvul,Rhoesc));
```

### Metazoa species set

List of species:

```tsv
# Species ID	Species name	Taxonomy	Data source
Nvec	Nematostella vectensis	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Actiniaria, Edwardsiidae	chromosomal assembly jaNemVect1.1, GCF_932526225.1; Combined JGI + Vienna annotation (https://doi.org/10.1016/j.cell.2018.05.019)
Exapal	Exaiptasia diaphana	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Actiniaria, Acuticulata, Aiptasiidae	NCBI GCF_001417965.1_Aiptasia_genome_1.1
Amil	Acropora millepora	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Acroporidae	NCBI assembly: GCA_013753865.1 (v2.1); annotation from Fuller 2020 (via https://drive.google.com/file/d/1ww7rbr6v8676Mul8OxoHkTp3fau0h27F/view?usp=sharing)
Spis	Stylophora pistillata	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Pocilloporidae	NCBI GCA_002571385.1 Stylophora pistillata v1 + reefgenomics (see Methods from Levy et al Cell 2021)
Ocupat	Oculina patagonica	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Oculinidae	This study
Ocuarb	Oculina arbuscula	Metazoa, Cnidaria, Anthozoa, Hexacorallia, Scleractinia, Oculinidae	NCBI GCA_964656845.1: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964656845.1/
Xesp	Xenia sp.	Metazoa, Cnidaria, Anthozoa, Octocorallia	Carnegie Coral & Marine Organisms database, https://cmo.carnegiescience.edu/data/
Dgig	Dendronephthya gigantea	Metazoa, Cnidaria, Anthozoa, Octocorallia	NCBI GCF_004324835.1_DenGig_1.0
Hvul	Hydra vulgaris	Metazoa, Cnidaria, Medusozoa, Hydrozoa, Anthoathecata, Hydridae	NIH https://research.nhgri.nih.gov/hydra/
Turdoh	Turritopsis dohrnii	Metazoa, Cnidaria, Medusozoa, Hydrozoa, Anthoathecata, Oceaniidae	Kazusa turritopsis.kazusa.or.jp, v2.0.1
Chem	Clytia hemisphaerica	Metazoa, Cnidaria, Medusozoa, Hydrozoa, Leptothecata, Clytiidae	MARIMBA Database http://marimba.obs-vlfr.fr/downloads/
Aaur	Aurelia aurita	Metazoa, Cnidaria, Medusozoa, Scyphozoa, Semaeostomeae	OIST ABSv1 (Atlantic strain, Baltic isolate) https://marinegenomics.oist.jp/aurelia_aurita/viewer/download?project_id=69
Rhoesc	Rhopilema esculentum	Metazoa, Cnidaria, Medusozoa, Scyphozoa, Rhizostomeae	NCBI GCA_013076305.1_ASM1307630v1 / annotation: https://github.com/nongwy/JellyfishGenomeData (scaffolds need to be mapped to NCBI names)
Hsap	Homo sapiens	Metazoa, Bilateria, Deuterostomia, Chordata, Vertebrata, Mammalia	Ensembl 102
Mmus	Mus musculus	Metazoa, Bilateria, Deuterostomia, Chordata, Vertebrata, Mammalia	Ensembl 102
Bralan	Branchiostoma lanceolatum	Metazoa, Bilateria, Deuterostomia, Chordata, Cephalochordata	NCBI GCA_927797965.1; v3
Spur	Strongylocentrotus purpuratus	Metazoa, Bilateria, Deuterostomia, Ambulacraria, Echinodermata, Echinoidea	NCBI GCA_000002235.4 Spur_5.0
Astrub	Asterias rubens	Metazoa, Bilateria, Deuterostomia, Ambulacraria, Echinodermata, Asteroidea	NCBI GCF_902459465.1_eAstRub1.3
Skow	Saccoglossus kowalevskii	Metazoa, Bilateria, Deuterostomia, Ambulacraria, Hemichordata	NCBI GCA_000003605.1 Skow_1.1
Dmel	Drosophila melanogaster	Metazoa, Bilateria, Protostomia, Ecdysozoa, Arthropoda, Insecta, Diptera	Ensembl Metazoa 49 BDGP6.28
Pricau	Priapulus caudatus	Metazoa, Bilateria, Protostomia, Ecdysozoa, Scalidophora, Priapulida	GCF_000485595.1 Priapulus caudatus-5.0.1
Owefus	Owenia fusiformis	Metazoa, Bilateria, Protostomia, Lophotrochozoa, Annelida	Liang et al. 2022, https://github.com/ChemaMD/OweniaGenome
Lgig	Lottia gigantea	Metazoa, Bilateria, Protostomia, Lophotrochozoa, Mollusca, Gastropoda	NCBI GCF_000327385.1_Helro1
Tadh	Trichoplax adhaerens H1	Metazoa, Placozoa	Ensembl Metazoa 49
Hhon	Hoilungia hongkongensis H13	Metazoa, Placozoa	Bitbucket: https://bitbucket.org/molpalmuc/hoilungia-genome/src/master/
```

Species tree:

```bash
((((((Nvec,Exapal)Actiniaria,(Amil,(Spis,Ocupat)OculPoci)Scleractinia)Hexacorallia,(Xesp,Dgig)Octocorallia)Anthozoa,(((Hvul,Turdoh)Anthoathecata,Chem)Hydrozoa,(Aaur,Rhoesc)Scyphozoa)Medusozoa)Cnidaria,((((Hsap,Mmus)Vertebrata,Bralan)Chordata,((Spur,Astrub)Echinodermata,Skow)Ambulacraria)Deuterostomia,((Dmel,Pricau)Ecdysozoa,(Owefus,Lgig)Lophotrochozoa)Protostomia)Bilateria)Planulozoa,(Tadh,Hhon)Placozoa)Parahoxozoa;
```
