
-- remove wrong domain annotations
update variant_v14 set pfam_name=null, pfam_description=null;
update variant_observation_v14 set pfam_name=null, pfam_description=null;
update subclonal_variant_observation_v14 set pfam_name=null, pfam_description=null;

-- set the right annotations on the variants table
select distinct on (d.name, d.gene_name) d.name as pfam_name, d.gene_name, g.start + d.start as start, g.start + d.end as end
    from gene_test as g, domain_test as d
    where g.name = d.gene_name;

update variant_v14 set pfam_name='bCoV_NAR', pfam_description='Non-structural protein NSP3, nucleic acid-binding (NAR) domain, betacoronavirus' where position >= 2188 and position <= 2285;
update variant_v14 set pfam_name='bCoV_NS6', pfam_description='Non-structural protein NS6, betacoronavirus' where position >= 27203 and position <= 27263;
update variant_v14 set pfam_name='bCoV_NS7A', pfam_description='Non-structural protein NS7A, coronavirus' where position >= 27410 and position <= 27515;
update variant_v14 set pfam_name='bCoV_NS7B', pfam_description='Non-structural protein 7b, SARS-like' where position >= 27757 and position <= 27798;
update variant_v14 set pfam_name='bCoV_NS8', pfam_description='Non-structural protein NS8, betacoronavirus' where position >= 27895 and position <= 28012;
update variant_v14 set pfam_name='bCoV_NSP1', pfam_description='Non-structural protein NSP1, betacoronavirus' where position >= 275 and position <= 409;
update variant_v14 set pfam_name='bCoV_NSP3_N', pfam_description='Non-structural protein NSP3, N-terminal, betacoronavirus' where position >= 1146 and position <= 1316;
update variant_v14 set pfam_name='bCoV_S1_N', pfam_description='Betacoronavirus-like spike glycoprotein S1, N-terminal' where position >= 21596 and position <= 21900;
update variant_v14 set pfam_name='bCoV_S1_RBD', pfam_description='Spike receptor binding domain, betacoronavirus' where position >= 21912 and position <= 22089;
update variant_v14 set pfam_name='bCoV_SUD_C', pfam_description='DPUP/SUD, C-terminal, betacoronavirus' where position >= 1764 and position <= 1827;
update variant_v14 set pfam_name='bCoV_SUD_M', pfam_description='Non-structural protein NSP3, single-stranded poly(A) binding domain, betacoronavirus' where position >= 1617 and position <= 1759;
update variant_v14 set pfam_name='bCoV_viroporin', pfam_description='Protein 3a, betacoronavirus' where position >= 25394 and position <= 25667;
update variant_v14 set pfam_name='CoV_E', pfam_description='Envelope small membrane protein, coronavirus' where position >= 26248 and position <= 26313;
update variant_v14 set pfam_name='CoV_M', pfam_description='M matrix/glycoprotein, coronavirus' where position >= 26539 and position <= 26739;
update variant_v14 set pfam_name='CoV_Methyltr_1', pfam_description='Non-structural protein NSP15, coronavirus' where position >= 6195 and position <= 6716;
update variant_v14 set pfam_name='CoV_Methyltr_2', pfam_description='Non-structural protein NSP16, coronavirus-like' where position >= 7066 and position <= 7361;
update variant_v14 set pfam_name='CoV_NSP10', pfam_description='RNA synthesis protein NSP10, coronavirus' where position >= 4528 and position <= 4650;
update variant_v14 set pfam_name='CoV_NSP15_C', pfam_description='Coronavirus replicase NSP15, uridylate-specific endoribonuclease' where position >= 6910 and position <= 7062;
update variant_v14 set pfam_name='CoV_NSP15_M', pfam_description='Coronavirus replicase NSP15, middle domain' where position >= 6789 and position <= 6885;
update variant_v14 set pfam_name='CoV_NSP15_N', pfam_description='Coronavirus replicase NSP15, N-terminal oligomerisation' where position >= 6719 and position <= 6779;
update variant_v14 set pfam_name='CoV_NSP2_C', pfam_description='Coronavirus replicase NSP2, C-terminal' where position >= 918 and position <= 1084;
update variant_v14 set pfam_name='CoV_NSP2_N', pfam_description='Coronavirus replicase NSP2, N-terminal' where position >= 449 and position <= 689;
update variant_v14 set pfam_name='CoV_NSP3_C', pfam_description='Coronavirus replicase NSP3, C-terminal' where position >= 2527 and position <= 3014;
update variant_v14 set pfam_name='CoV_NSP4_C', pfam_description='Non-structural protein NSP4, C-terminal, coronavirus' where position >= 3432 and position <= 3527;
update variant_v14 set pfam_name='CoV_NSP4_N', pfam_description='Coronavirus replicase NSP4, N-terminal' where position >= 3054 and position <= 3407;
update variant_v14 set pfam_name='CoV_NSP6', pfam_description='Coronavirus replicase NSP6' where position >= 3864 and position <= 4125;
update variant_v14 set pfam_name='CoV_NSP7', pfam_description='Non-structural protein NSP7, coronavirus' where position >= 4126 and position <= 4208;
update variant_v14 set pfam_name='CoV_NSP8', pfam_description='Non-structural protein NSP8, coronavirus-like' where position >= 4209 and position <= 4405;
update variant_v14 set pfam_name='CoV_NSP9', pfam_description='Non-structural protein NSP9, coronavirus' where position >= 4407 and position <= 4519;
update variant_v14 set pfam_name='CoV_nucleocap', pfam_description='Nucleocapsid protein, coronavirus' where position >= 28316 and position <= 28656;
update variant_v14 set pfam_name='CoV_peptidase', pfam_description='Peptidase C16, coronavirus' where position >= 1830 and position <= 2148;
update variant_v14 set pfam_name='CoV_RPol_N', pfam_description='RNA polymerase, N-terminal, coronaviral' where position >= 4673 and position <= 5024;
update variant_v14 set pfam_name='CoV_S1_C', pfam_description='Coronavirus spike glycoprotein S1, C-terminal' where position >= 22099 and position <= 22155;
update variant_v14 set pfam_name='CoV_S2', pfam_description='Spike glycoprotein S2, coronavirus' where position >= 22274 and position <= 22795;
update variant_v14 set pfam_name='Macro', pfam_description='Macro domain' where position >= 1324 and position <= 1430;
update variant_v14 set pfam_name='Peptidase_C30', pfam_description='Peptidase C30, coronavirus' where position >= 3558 and position <= 3848;
update variant_v14 set pfam_name='RdRP_1', pfam_description='RNA-directed RNA polymerase,  C-terminal domain' where position >= 5058 and position <= 5546;
update variant_v14 set pfam_name='Viral_helicase1', pfam_description='(+) RNA virus helicase core domain' where position >= 5938 and position <= 6162;

-- amend other tables from the variants table
update variant_observation_v14
    set pfam_name=v.pfam_name, pfam_description=v.pfam_description
    from variant_v14 as v where v.variant_id = variant_observation_v14.variant_id;
update subclonal_variant_observation_v14
    set pfam_name=v.pfam_name, pfam_description=v.pfam_description
    from variant_v14 as v where v.variant_id = subclonal_variant_observation_v14.variant_id;