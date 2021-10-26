
-- remove wrong domain annotations
update variant_v19 set pfam_name=null, pfam_description=null;

-- gets the right annotations on the variants table
select distinct on (d.name, d.gene_name) d.name as pfam_name, d.gene_name, g.start + (d.start * 3) as start, g.start + (d.end * 3) as end
    from gene_test as g, domain_test as d
    where g.name = d.gene_name;

-- sets the right annotations in variants
update variant_v19 set pfam_name='bCoV_NAR', pfam_description='Non-structural protein NSP3, nucleic acid-binding (NAR) domain, betacoronavirus' where position >= 6032 and position <= 6323;
update variant_v19 set pfam_name='bCoV_NS6', pfam_description='Non-structural protein NS6, betacoronavirus' where position >= 27205 and position <= 27385;
update variant_v19 set pfam_name='bCoV_NS7A', pfam_description='Non-structural protein NS7A, coronavirus' where position >= 27442 and position <= 27757;
update variant_v19 set pfam_name='bCoV_NS7B', pfam_description='Non-structural protein 7b, SARS-like' where position >= 27759 and position <= 27882;
update variant_v19 set pfam_name='bCoV_NS8', pfam_description='Non-structural protein NS8, betacoronavirus' where position >= 27897 and position <= 28248;
update variant_v19 set pfam_name='bCoV_NSP1', pfam_description='Non-structural protein NSP1, betacoronavirus' where position >= 293 and position <= 695;
update variant_v19 set pfam_name='bCoV_NSP3_N', pfam_description='Non-structural protein NSP3, N-terminal, betacoronavirus' where position >= 2906 and position <= 3416;
update variant_v19 set pfam_name='bCoV_S1_N', pfam_description='Betacoronavirus-like spike glycoprotein S1, N-terminal' where position >= 21662 and position <= 22574;
update variant_v19 set pfam_name='bCoV_S1_RBD', pfam_description='Spike receptor binding domain, betacoronavirus' where position >= 22610 and position <= 23141;
update variant_v19 set pfam_name='bCoV_SUD_C', pfam_description='DPUP/SUD, C-terminal, betacoronavirus' where position >= 4760 and position <= 4949;
update variant_v19 set pfam_name='bCoV_SUD_M', pfam_description='Non-structural protein NSP3, single-stranded poly(A) binding domain, betacoronavirus' where position >= 4319 and position <= 4745;
update variant_v19 set pfam_name='bCoV_viroporin', pfam_description='Protein 3a, betacoronavirus' where position >= 25396 and position <= 26215;
update variant_v19 set pfam_name='CoV_E', pfam_description='Envelope small membrane protein, coronavirus' where position >= 26254 and position <= 26449;
update variant_v19 set pfam_name='CoV_M', pfam_description='M matrix/glycoprotein, coronavirus' where position >= 26571 and position <= 27171;
update variant_v19 set pfam_name='CoV_Methyltr_1', pfam_description='Non-structural protein NSP15, coronavirus' where position >= 18053 and position <= 19616;
update variant_v19 set pfam_name='CoV_Methyltr_2', pfam_description='Non-structural protein NSP16, coronavirus-like' where position >= 20666 and position <= 21551;
update variant_v19 set pfam_name='CoV_NSP10', pfam_description='RNA synthesis protein NSP10, coronavirus' where position >= 13052 and position <= 13418;
update variant_v19 set pfam_name='CoV_NSP15_C', pfam_description='Coronavirus replicase NSP15, uridylate-specific endoribonuclease' where position >= 20198 and position <= 20654;
update variant_v19 set pfam_name='CoV_NSP15_M', pfam_description='Coronavirus replicase NSP15, middle domain' where position >= 19835 and position <= 20123;
update variant_v19 set pfam_name='CoV_NSP15_N', pfam_description='Coronavirus replicase NSP15, N-terminal oligomerisation' where position >= 19625 and position <= 19805;
update variant_v19 set pfam_name='CoV_NSP2_C', pfam_description='Coronavirus replicase NSP2, C-terminal' where position >= 2222 and position <= 2720;
update variant_v19 set pfam_name='CoV_NSP2_N', pfam_description='Coronavirus replicase NSP2, N-terminal' where position >= 815 and position <= 1535;
update variant_v19 set pfam_name='CoV_NSP3_C', pfam_description='Coronavirus replicase NSP3, C-terminal' where position >= 7049 and position <= 8510;
update variant_v19 set pfam_name='CoV_NSP4_C', pfam_description='Non-structural protein NSP4, C-terminal, coronavirus' where position >= 9764 and position <= 10049;
update variant_v19 set pfam_name='CoV_NSP4_N', pfam_description='Coronavirus replicase NSP4, N-terminal' where position >= 8630 and position <= 9689;
update variant_v19 set pfam_name='CoV_NSP6', pfam_description='Coronavirus replicase NSP6' where position >= 11060 and position <= 11843;
update variant_v19 set pfam_name='CoV_NSP7', pfam_description='Non-structural protein NSP7, coronavirus' where position >= 11846 and position <= 12092;
update variant_v19 set pfam_name='CoV_NSP8', pfam_description='Non-structural protein NSP8, coronavirus-like' where position >= 12095 and position <= 12683;
update variant_v19 set pfam_name='CoV_NSP9', pfam_description='Non-structural protein NSP9, coronavirus' where position >= 12689 and position <= 13025;
update variant_v19 set pfam_name='CoV_nucleocap', pfam_description='Nucleocapsid protein, coronavirus' where position >= 28400 and position <= 29420;
update variant_v19 set pfam_name='CoV_peptidase', pfam_description='Peptidase C16, coronavirus' where position >= 4958 and position <= 5912;
update variant_v19 set pfam_name='CoV_RPol_N', pfam_description='RNA polymerase, N-terminal, coronaviral' where position >= 13487 and position <= 14540;
update variant_v19 set pfam_name='CoV_S1_C', pfam_description='Coronavirus spike glycoprotein S1, C-terminal' where position >= 23171 and position <= 23339;
update variant_v19 set pfam_name='CoV_S2', pfam_description='Spike glycoprotein S2, coronavirus' where position >= 23696 and position <= 25259;
update variant_v19 set pfam_name='Macro', pfam_description='Macro domain' where position >= 3440 and position <= 3758;
update variant_v19 set pfam_name='Peptidase_C30', pfam_description='Peptidase C30, coronavirus' where position >= 10142 and position <= 11012;
update variant_v19 set pfam_name='RdRP_1', pfam_description='RNA-directed RNA polymerase,  C-terminal domain' where position >= 14642 and position <= 16106;
update variant_v19 set pfam_name='Viral_helicase1', pfam_description='(+) RNA virus helicase core domain' where position >= 17282 and position <= 17954;


-- sets the right annotations in variant_observations
ALTER TABLE variant_observation_v19 RENAME COLUMN pfam_name TO old_pfam_name;
ALTER TABLE variant_observation_v19 RENAME COLUMN pfam_description TO old_pfam_description;
ALTER TABLE variant_observation_v19 ADD column pfam_name varchar;
ALTER TABLE variant_observation_v19 ADD column pfam_description varchar;

update variant_observation_v19 set pfam_name='bCoV_NAR', pfam_description='Non-structural protein NSP3, nucleic acid-binding (NAR) domain, betacoronavirus' where position >= 6032 and position <= 6323;
update variant_observation_v19 set pfam_name='bCoV_NS6', pfam_description='Non-structural protein NS6, betacoronavirus' where position >= 27205 and position <= 27385;
update variant_observation_v19 set pfam_name='bCoV_NS7A', pfam_description='Non-structural protein NS7A, coronavirus' where position >= 27442 and position <= 27757;
update variant_observation_v19 set pfam_name='bCoV_NS7B', pfam_description='Non-structural protein 7b, SARS-like' where position >= 27759 and position <= 27882;
update variant_observation_v19 set pfam_name='bCoV_NS8', pfam_description='Non-structural protein NS8, betacoronavirus' where position >= 27897 and position <= 28248;
update variant_observation_v19 set pfam_name='bCoV_NSP1', pfam_description='Non-structural protein NSP1, betacoronavirus' where position >= 293 and position <= 695;
update variant_observation_v19 set pfam_name='bCoV_NSP3_N', pfam_description='Non-structural protein NSP3, N-terminal, betacoronavirus' where position >= 2906 and position <= 3416;
update variant_observation_v19 set pfam_name='bCoV_S1_N', pfam_description='Betacoronavirus-like spike glycoprotein S1, N-terminal' where position >= 21662 and position <= 22574;
update variant_observation_v19 set pfam_name='bCoV_S1_RBD', pfam_description='Spike receptor binding domain, betacoronavirus' where position >= 22610 and position <= 23141;
update variant_observation_v19 set pfam_name='bCoV_SUD_C', pfam_description='DPUP/SUD, C-terminal, betacoronavirus' where position >= 4760 and position <= 4949;
update variant_observation_v19 set pfam_name='bCoV_SUD_M', pfam_description='Non-structural protein NSP3, single-stranded poly(A) binding domain, betacoronavirus' where position >= 4319 and position <= 4745;
update variant_observation_v19 set pfam_name='bCoV_viroporin', pfam_description='Protein 3a, betacoronavirus' where position >= 25396 and position <= 26215;
update variant_observation_v19 set pfam_name='CoV_E', pfam_description='Envelope small membrane protein, coronavirus' where position >= 26254 and position <= 26449;
update variant_observation_v19 set pfam_name='CoV_M', pfam_description='M matrix/glycoprotein, coronavirus' where position >= 26571 and position <= 27171;
update variant_observation_v19 set pfam_name='CoV_Methyltr_1', pfam_description='Non-structural protein NSP15, coronavirus' where position >= 18053 and position <= 19616;
update variant_observation_v19 set pfam_name='CoV_Methyltr_2', pfam_description='Non-structural protein NSP16, coronavirus-like' where position >= 20666 and position <= 21551;
update variant_observation_v19 set pfam_name='CoV_NSP10', pfam_description='RNA synthesis protein NSP10, coronavirus' where position >= 13052 and position <= 13418;
update variant_observation_v19 set pfam_name='CoV_NSP15_C', pfam_description='Coronavirus replicase NSP15, uridylate-specific endoribonuclease' where position >= 20198 and position <= 20654;
update variant_observation_v19 set pfam_name='CoV_NSP15_M', pfam_description='Coronavirus replicase NSP15, middle domain' where position >= 19835 and position <= 20123;
update variant_observation_v19 set pfam_name='CoV_NSP15_N', pfam_description='Coronavirus replicase NSP15, N-terminal oligomerisation' where position >= 19625 and position <= 19805;
update variant_observation_v19 set pfam_name='CoV_NSP2_C', pfam_description='Coronavirus replicase NSP2, C-terminal' where position >= 2222 and position <= 2720;
update variant_observation_v19 set pfam_name='CoV_NSP2_N', pfam_description='Coronavirus replicase NSP2, N-terminal' where position >= 815 and position <= 1535;
update variant_observation_v19 set pfam_name='CoV_NSP3_C', pfam_description='Coronavirus replicase NSP3, C-terminal' where position >= 7049 and position <= 8510;
update variant_observation_v19 set pfam_name='CoV_NSP4_C', pfam_description='Non-structural protein NSP4, C-terminal, coronavirus' where position >= 9764 and position <= 10049;
update variant_observation_v19 set pfam_name='CoV_NSP4_N', pfam_description='Coronavirus replicase NSP4, N-terminal' where position >= 8630 and position <= 9689;
update variant_observation_v19 set pfam_name='CoV_NSP6', pfam_description='Coronavirus replicase NSP6' where position >= 11060 and position <= 11843;
update variant_observation_v19 set pfam_name='CoV_NSP7', pfam_description='Non-structural protein NSP7, coronavirus' where position >= 11846 and position <= 12092;
update variant_observation_v19 set pfam_name='CoV_NSP8', pfam_description='Non-structural protein NSP8, coronavirus-like' where position >= 12095 and position <= 12683;
update variant_observation_v19 set pfam_name='CoV_NSP9', pfam_description='Non-structural protein NSP9, coronavirus' where position >= 12689 and position <= 13025;
update variant_observation_v19 set pfam_name='CoV_nucleocap', pfam_description='Nucleocapsid protein, coronavirus' where position >= 28400 and position <= 29420;
update variant_observation_v19 set pfam_name='CoV_peptidase', pfam_description='Peptidase C16, coronavirus' where position >= 4958 and position <= 5912;
update variant_observation_v19 set pfam_name='CoV_RPol_N', pfam_description='RNA polymerase, N-terminal, coronaviral' where position >= 13487 and position <= 14540;
update variant_observation_v19 set pfam_name='CoV_S1_C', pfam_description='Coronavirus spike glycoprotein S1, C-terminal' where position >= 23171 and position <= 23339;
update variant_observation_v19 set pfam_name='CoV_S2', pfam_description='Spike glycoprotein S2, coronavirus' where position >= 23696 and position <= 25259;
update variant_observation_v19 set pfam_name='Macro', pfam_description='Macro domain' where position >= 3440 and position <= 3758;
update variant_observation_v19 set pfam_name='Peptidase_C30', pfam_description='Peptidase C30, coronavirus' where position >= 10142 and position <= 11012;
update variant_observation_v19 set pfam_name='RdRP_1', pfam_description='RNA-directed RNA polymerase,  C-terminal domain' where position >= 14642 and position <= 16106;
update variant_observation_v19 set pfam_name='Viral_helicase1', pfam_description='(+) RNA virus helicase core domain' where position >= 17282 and position <= 17954;


-- sets the right annotations in subclonal_variant_observations
ALTER TABLE subclonal_variant_observation_v19 RENAME COLUMN pfam_name TO old_pfam_name;
ALTER TABLE subclonal_variant_observation_v19 RENAME COLUMN pfam_description TO old_pfam_description;
ALTER TABLE subclonal_variant_observation_v19 ADD column pfam_name varchar;
ALTER TABLE subclonal_variant_observation_v19 ADD column pfam_description varchar;

update subclonal_variant_observation_v19 set pfam_name='bCoV_NAR', pfam_description='Non-structural protein NSP3, nucleic acid-binding (NAR) domain, betacoronavirus' where position >= 6032 and position <= 6323;
update subclonal_variant_observation_v19 set pfam_name='bCoV_NS6', pfam_description='Non-structural protein NS6, betacoronavirus' where position >= 27205 and position <= 27385;
update subclonal_variant_observation_v19 set pfam_name='bCoV_NS7A', pfam_description='Non-structural protein NS7A, coronavirus' where position >= 27442 and position <= 27757;
update subclonal_variant_observation_v19 set pfam_name='bCoV_NS7B', pfam_description='Non-structural protein 7b, SARS-like' where position >= 27759 and position <= 27882;
update subclonal_variant_observation_v19 set pfam_name='bCoV_NS8', pfam_description='Non-structural protein NS8, betacoronavirus' where position >= 27897 and position <= 28248;
update subclonal_variant_observation_v19 set pfam_name='bCoV_NSP1', pfam_description='Non-structural protein NSP1, betacoronavirus' where position >= 293 and position <= 695;
update subclonal_variant_observation_v19 set pfam_name='bCoV_NSP3_N', pfam_description='Non-structural protein NSP3, N-terminal, betacoronavirus' where position >= 2906 and position <= 3416;
update subclonal_variant_observation_v19 set pfam_name='bCoV_S1_N', pfam_description='Betacoronavirus-like spike glycoprotein S1, N-terminal' where position >= 21662 and position <= 22574;
update subclonal_variant_observation_v19 set pfam_name='bCoV_S1_RBD', pfam_description='Spike receptor binding domain, betacoronavirus' where position >= 22610 and position <= 23141;
update subclonal_variant_observation_v19 set pfam_name='bCoV_SUD_C', pfam_description='DPUP/SUD, C-terminal, betacoronavirus' where position >= 4760 and position <= 4949;
update subclonal_variant_observation_v19 set pfam_name='bCoV_SUD_M', pfam_description='Non-structural protein NSP3, single-stranded poly(A) binding domain, betacoronavirus' where position >= 4319 and position <= 4745;
update subclonal_variant_observation_v19 set pfam_name='bCoV_viroporin', pfam_description='Protein 3a, betacoronavirus' where position >= 25396 and position <= 26215;
update subclonal_variant_observation_v19 set pfam_name='CoV_E', pfam_description='Envelope small membrane protein, coronavirus' where position >= 26254 and position <= 26449;
update subclonal_variant_observation_v19 set pfam_name='CoV_M', pfam_description='M matrix/glycoprotein, coronavirus' where position >= 26571 and position <= 27171;
update subclonal_variant_observation_v19 set pfam_name='CoV_Methyltr_1', pfam_description='Non-structural protein NSP15, coronavirus' where position >= 18053 and position <= 19616;
update subclonal_variant_observation_v19 set pfam_name='CoV_Methyltr_2', pfam_description='Non-structural protein NSP16, coronavirus-like' where position >= 20666 and position <= 21551;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP10', pfam_description='RNA synthesis protein NSP10, coronavirus' where position >= 13052 and position <= 13418;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP15_C', pfam_description='Coronavirus replicase NSP15, uridylate-specific endoribonuclease' where position >= 20198 and position <= 20654;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP15_M', pfam_description='Coronavirus replicase NSP15, middle domain' where position >= 19835 and position <= 20123;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP15_N', pfam_description='Coronavirus replicase NSP15, N-terminal oligomerisation' where position >= 19625 and position <= 19805;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP2_C', pfam_description='Coronavirus replicase NSP2, C-terminal' where position >= 2222 and position <= 2720;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP2_N', pfam_description='Coronavirus replicase NSP2, N-terminal' where position >= 815 and position <= 1535;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP3_C', pfam_description='Coronavirus replicase NSP3, C-terminal' where position >= 7049 and position <= 8510;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP4_C', pfam_description='Non-structural protein NSP4, C-terminal, coronavirus' where position >= 9764 and position <= 10049;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP4_N', pfam_description='Coronavirus replicase NSP4, N-terminal' where position >= 8630 and position <= 9689;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP6', pfam_description='Coronavirus replicase NSP6' where position >= 11060 and position <= 11843;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP7', pfam_description='Non-structural protein NSP7, coronavirus' where position >= 11846 and position <= 12092;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP8', pfam_description='Non-structural protein NSP8, coronavirus-like' where position >= 12095 and position <= 12683;
update subclonal_variant_observation_v19 set pfam_name='CoV_NSP9', pfam_description='Non-structural protein NSP9, coronavirus' where position >= 12689 and position <= 13025;
update subclonal_variant_observation_v19 set pfam_name='CoV_nucleocap', pfam_description='Nucleocapsid protein, coronavirus' where position >= 28400 and position <= 29420;
update subclonal_variant_observation_v19 set pfam_name='CoV_peptidase', pfam_description='Peptidase C16, coronavirus' where position >= 4958 and position <= 5912;
update subclonal_variant_observation_v19 set pfam_name='CoV_RPol_N', pfam_description='RNA polymerase, N-terminal, coronaviral' where position >= 13487 and position <= 14540;
update subclonal_variant_observation_v19 set pfam_name='CoV_S1_C', pfam_description='Coronavirus spike glycoprotein S1, C-terminal' where position >= 23171 and position <= 23339;
update subclonal_variant_observation_v19 set pfam_name='CoV_S2', pfam_description='Spike glycoprotein S2, coronavirus' where position >= 23696 and position <= 25259;
update subclonal_variant_observation_v19 set pfam_name='Macro', pfam_description='Macro domain' where position >= 3440 and position <= 3758;
update subclonal_variant_observation_v19 set pfam_name='Peptidase_C30', pfam_description='Peptidase C30, coronavirus' where position >= 10142 and position <= 11012;
update subclonal_variant_observation_v19 set pfam_name='RdRP_1', pfam_description='RNA-directed RNA polymerase,  C-terminal domain' where position >= 14642 and position <= 16106;
update subclonal_variant_observation_v19 set pfam_name='Viral_helicase1', pfam_description='(+) RNA virus helicase core domain' where position >= 17282 and position <= 17954;



-- THIS IS TOO SLOW...
-- amend other tables from the variants table
alter table variant_observation_v19 set (FILLFACTOR=70);
VACUUM FULL variant_observation_v19;
REINDEX TABLE variant_observation_v19;
update variant_observation_v19
    set pfam_name=v.pfam_name, pfam_description=v.pfam_description
    from variant_v19 as v where v.variant_id = variant_observation_v19.variant_id;

alter table subclonal_variant_observation_v19 set (FILLFACTOR=70);
VACUUM FULL subclonal_variant_observation_v19;
REINDEX TABLE subclonal_variant_observation_v19;
update subclonal_variant_observation_v19
    set pfam_name=v.pfam_name, pfam_description=v.pfam_description
    from variant_v19 as v where v.variant_id = subclonal_variant_observation_v19.variant_id;


-- TESTING FASTER THINGS LOCALLY
CREATE TABLE variant_observation_ena_v19 AS
  SELECT *
  FROM variant_observation_v19
  WHERE source = 'ENA';

