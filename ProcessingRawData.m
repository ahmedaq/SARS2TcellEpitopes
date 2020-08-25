clear all;close all;clc

% run startup.m

%%

t_table = readtable('Raw Data/SARS-CoV-2-T-cell-epitopes-v2.xlsx');

t_table.Protein(find(strcmp(t_table.Protein,'orf1a')))={'ORF1a'};
t_table.Protein(find(strcmp(t_table.Protein,'orf1b')))={'ORF1b'};

unique_proteins = unique(t_table.Protein);


%% Proteins

peptides_all = t_table.Peptide;

[peptides_unique,bb,cc] = unique(peptides_all);

for kk = 1:length(bb)
    peptides_unique_indices{kk} = find(cc == cc(bb(kk)));
    %     no_seqs_in_each_country(kk) = length(find(c == c(b(kk))));
    peptides_unique_count(kk) = length(peptides_unique_indices{kk});
end


unique_proteins = unique(t_table.Protein);

% for kk = 1:length(peptides_unique_indices)
%     peptides_unique_indices2(kk) = peptides_unique_indices{kk}(1);
% end
% 
% for kk = 1:length(unique_proteins)
%     no_epitopes_protein(kk) = sum(strcmp(t_table.Protein(peptides_unique_indices2),unique_proteins{kk}));
% end


%%

t_table.No_OfSubjectsResponded = str2double(t_table.No_OfSubjectsResponded);
t_table.No_OfTotalSubjects = str2double(t_table.No_OfTotalSubjects);
t_table.Start = str2double(t_table.Start);
t_table.Stop = str2double(t_table.Stop);

%% Calculating RF [CI]

num_responded = zeros(1,length(peptides_unique_indices));
num_total = zeros(1,length(peptides_unique_indices));

rf = zeros(length(peptides_unique_indices),1);
ci = zeros(length(peptides_unique_indices),2);

for kk = 1:length(peptides_unique_indices)
    
    if length(peptides_unique_indices{kk})>1 %if multiple entries for one epitope
        
        if sum(isnan(t_table.No_OfSubjectsResponded(peptides_unique_indices{kk}))) == length(peptides_unique_indices{kk}) %all nan
            num_responded(kk) = nan;
            num_total(kk) = nan;
            
        elseif sum(isnan(t_table.No_OfSubjectsResponded(peptides_unique_indices{kk}))) <= length(peptides_unique_indices{kk}) %at least one value
            indices_data = find(~isnan(t_table.No_OfSubjectsResponded(peptides_unique_indices{kk})));
            
            for mm = 1:length(indices_data)
                num_responded(kk) = num_responded(kk) + (t_table.No_OfSubjectsResponded(peptides_unique_indices{kk}(indices_data(mm))));
                num_total(kk) = num_total(kk) + (t_table.No_OfTotalSubjects(peptides_unique_indices{kk}(indices_data(mm))));
            end
            
        else %all values
            for mm = 1:length(peptides_unique_indices{kk})
                
                num_responded(kk) = num_responded(kk) + (t_table.No_OfSubjectsResponded(peptides_unique_indices{kk}(mm)));
                num_total(kk) = num_total(kk) + (t_table.No_OfTotalSubjects(peptides_unique_indices{kk}(mm)));
                
            end
            
        end
        
    else % one entry for an epitope
        
        num_responded(kk) = (t_table.No_OfSubjectsResponded(peptides_unique_indices{kk}(1)));
        num_total(kk) = num_total(kk) + (t_table.No_OfTotalSubjects(peptides_unique_indices{kk}(1)));
        
    end
    
    if sum(isnan(num_responded(kk)))
        
        rf(kk) = nan;
        ci(kk,1:2) = nan;
        
    else
        
        [rf(kk),ci(kk,:)] = binofit(num_responded(kk),num_total(kk));
    end
    
end


%% putting rf in the t_table

%make a table of unique epitopes of SARS-CoV-2 and add rf to it


% t_table_rf = array2table(zeros(0,size(t_table,2)));
% t_table_rf.Properties = t_table.Properties;

t_table_rf = table('Size',[length(peptides_unique) size(t_table,2)],...
    'VariableTypes',{'cell','double','double','cell','cell','cell','double','double','cell','cell','cell','cell','cell','double'},...
    'VariableNames',t_table.Properties.VariableNames);


for kk = 1:length(peptides_unique)
    
    t_table_rf.Protein{kk} = t_table.Protein{peptides_unique_indices{kk}(1)};
    t_table_rf.Start(kk) = t_table.Start(peptides_unique_indices{kk}(1));
    t_table_rf.Stop(kk) = t_table.Stop(peptides_unique_indices{kk}(1));
    t_table_rf.Peptide{kk} = t_table.Peptide(peptides_unique_indices{kk}(1));
    t_table_rf.Response_Y_N_{kk} = t_table.Response_Y_N_{peptides_unique_indices{kk}(1)};
    t_table_rf.No_OfSubjectsResponded(kk) = num_responded(kk);
    t_table_rf.No_OfTotalSubjects(kk) = num_total(kk);
    t_table_rf.ResponseFreq(kk) = rf(kk);
    t_table_rf.ResponseFreq_CI_min(kk) = ci(kk,1);
    t_table_rf.ResponseFreq_CI_max(kk) = ci(kk,2);
    t_table_rf.Conservation(kk) = t_table.Conservation(peptides_unique_indices{kk}(1));
    
    %         if length(peptides_unique_indices)>1
    
    
    
    for mm = 1:length(peptides_unique_indices)
        
        unique_host = unique(t_table.Host(peptides_unique_indices{kk}));
        for nn = 1:length(unique_host)
            if nn == 1
                t_table_rf.Host{kk} = unique_host{nn};
            else
                t_table_rf.Host{kk} = ...
                    sprintf('%s/%s',t_table_rf.Host{kk},unique_host{nn});
            end
        end
        
        
        unique_HLAs = unique(t_table.HLAAlleleInformation__Separated_(peptides_unique_indices{kk}));
        for nn = 1:length(unique_HLAs)
            if nn == 1
                t_table_rf.HLAAlleleInformation__Separated_{kk} = unique_HLAs{nn};
                t_table_rf.HLAClassInformation_I_II_{kk} = t_table.HLAClassInformation_I_II_{peptides_unique_indices{kk}(nn)};
            else
                t_table_rf.HLAAlleleInformation__Separated_{kk} = ...
                    sprintf('%s/%s',t_table_rf.HLAAlleleInformation__Separated_{kk},unique_HLAs{nn});
                t_table_rf.HLAClassInformation_I_II_{kk} = ...
                    sprintf('%s/%s',t_table_rf.HLAClassInformation_I_II_{kk},t_table.HLAClassInformation_I_II_{peptides_unique_indices{kk}(nn)});
            end
        end
        
        unique_doi = unique(t_table.PublicationID_doi_(peptides_unique_indices{kk}));
        for nn = 1:length(unique_doi)
            if nn == 1
                t_table_rf.PublicationID_doi_{kk} = unique_doi{nn};
                t_table_rf.PublicationRef_{kk} = t_table.PublicationRef_{peptides_unique_indices{kk}(nn)};
                t_table_rf.Comments_Notes{kk} = t_table.Comments_Notes{peptides_unique_indices{kk}(nn)};
                
            else
                t_table_rf.PublicationID_doi_{kk} = ...
                    sprintf('%s,%s',t_table_rf.PublicationID_doi_{kk},unique_doi{nn});
                t_table_rf.PublicationRef_{kk} = ...
                    sprintf('%s/%s',t_table_rf.PublicationRef_{kk},t_table.PublicationRef_{peptides_unique_indices{kk}(nn)});
                t_table_rf.Comments_Notes{kk} = ...
                    sprintf('%s/%s',t_table_rf.Comments_Notes{kk},t_table.Comments_Notes{peptides_unique_indices{kk}(nn)});
                
            end
        end
        
    end
        
   kk     
end

writetable(t_table_rf,'Processed Data/SARS-CoV-2-T-cell-epitopes-RespFreq.xlsx')



%% Length of epitopes


for kk = 1:length(peptides_unique)
    length_peptides(kk) = length(peptides_unique{kk});%length(t_table.Start(kk):t_table.Stop(kk));
end

%% HLA allele information of epitopes

no_HLAs = zeros(1,length(peptides_unique_indices));
hla_info_peptide = cell(1,length(peptides_unique_indices));

for kk = 1:length(peptides_unique_indices)
    
%     hla_info_peptide{kk} = [];
    
    if length(peptides_unique_indices{kk})>1
       
        unique_HLAs = unique(t_table.HLAAlleleInformation__Separated_(peptides_unique_indices{kk}));
        for nn = 1:length(unique_HLAs)
            if nn == 1
                hla_info_peptide{kk} = unique_HLAs{nn};
            else
                hla_info_peptide{kk} = ...
                    sprintf('%s/%s',t_table_rf.HLAAlleleInformation__Separated_{kk},unique_HLAs{nn});
            end
        end
        
    else
        
        hla_info_peptide{kk} = t_table.HLAAlleleInformation__Separated_{peptides_unique_indices{kk}};
        
    end
    
    
    if isempty(hla_info_peptide{kk})
        no_HLAs(kk) = 0;
    else
        no_HLAs(kk) = length(find(hla_info_peptide{kk} == '/')) + 1;
    end
end



%% Loading SARS-CoV epitope data downloaded from COVIDep

t_table3 = readtable('Raw Data/AllSARSTepitopes.csv');
% t_table3 = readtable('Raw Data/COVIDepALLTcell.csv');


for kk = 1:size(t_table3,1)
    
    alleles = t_table3.MHCAlleleNames(kk);
%     t_table2.MHCAlleleNames(kk) = alleles;
    if strcmp(alleles,'-N/A-') || strcmp(alleles,'HLA class II') || strcmp(alleles,'HLA class I')
%         t_table2.NoMHCs(kk) = 0;
        t_table3.NumMHCs(kk) = 0;
    else
%         t_table2.NoMHCs(kk) = sum(alleles{:}=='/') + 1;
        t_table3.NumMHCs(kk) = sum(alleles{:}=='/') + 1;
    end
end


%% Finding epitopes completely overlapping the peptides in t_table_rf

t_table_mapping = table('Size',[20 20],...
    'VariableTypes',{'cell','cell','cell','cell','cell','cell','double',...
    'cell','cell','cell','cell','double','double','double','double','double'...
    'cell','cell','cell','cell'},...
    'VariableNames',{'IEDB ID','Protein','SARS Epitope','SARS T Cell Assay','SARS HLA Class','SARS HLA Alleles','SARS Number of HLAs',...
    'SARS2 Epitope/Peptide','SARS2 Host','SARS2 HLA Alleles','SARS2 HLA Class','SARS2 Number of Responses','SARS2 Total Subjects','RF','RF CI min','RF CI max'...
    'Epitope (refined)','doi','Reference','Comments'});


nn = 1;
for kk = 1:size(t_table_rf,1)
    
    nnn = nn;
    
    for mm = 1:size(t_table3,1)
        
        [aa,bb,cc] = swalign(t_table_rf.Peptide{kk}{:},t_table3.Epitope{mm});
        
        if isempty(t_table_rf.HLAAlleleInformation__Separated_{kk})
        
            perc_overlap = sum(bb(2,:)=='|')/length(t_table3.Epitope{mm})*100; %allows fragmets of sars2 to match sars epitope
            
        else
            
            perc_overlap = sum(bb(2,:)=='|')*2/(length(t_table_rf.Peptide{kk}{:})+length(t_table3.Epitope{mm}))*100; %exact match only (of same length)
            
        end
        
        if perc_overlap==100 %when 100% overlap
            
            if ~isempty(t_table_rf.HLAClassInformation_I_II_{kk}) %if HLA class info is NOT empty
                
                if strcmp(t_table_rf.HLAClassInformation_I_II_{kk},t_table3.MHCAlleleClass{mm}) %if HLA class matches
                    
                    t_table_mapping.("IEDB ID"){nn} = t_table3.IEDB(mm);
                    t_table_mapping.Protein{nn} = t_table_rf.Protein{kk};
                    t_table_mapping.("SARS Epitope"){nn} = t_table3.Epitope{mm};
%                     t_table_mapping.("SARS T Cell Assay"){nn} = t_table3.TCellAssay{mm};
                    t_table_mapping.("SARS HLA Class"){nn} = t_table3.MHCAlleleClass{mm};
                    t_table_mapping.("SARS HLA Alleles"){nn} = t_table3.MHCAlleleNames{mm};
                    t_table_mapping.("SARS Number of HLAs")(nn) = t_table3.NumMHCs(mm);
                    t_table_mapping.("SARS2 Epitope/Peptide"){nn} = t_table_rf.Peptide{kk}{:};
                    t_table_mapping.("SARS2 Host"){nn} = t_table_rf.Host{kk};
                    t_table_mapping.("SARS2 HLA Alleles"){nn} = t_table_rf.HLAAlleleInformation__Separated_{kk};
                    t_table_mapping.("SARS2 HLA Class"){nn} = t_table_rf.HLAClassInformation_I_II_{kk};
                    t_table_mapping.("SARS2 Number of Responses")(nn) = t_table_rf.No_OfSubjectsResponded(kk);
                    t_table_mapping.("SARS2 Total Subjects")(nn) = t_table_rf.No_OfSubjectsResponded(kk);
                    t_table_mapping.("RF")(nn) = t_table_rf.ResponseFreq(kk);
                    t_table_mapping.("RF CI min")(nn) = t_table_rf.ResponseFreq_CI_min(kk);
                    t_table_mapping.("RF CI max")(nn) = t_table_rf.ResponseFreq_CI_max(kk);
                    
                    t_table_mapping.("Epitope (refined)"){nn} = t_table3.Epitope{mm};
                    t_table_mapping.("doi"){nn} = t_table_rf.PublicationID_doi_{kk};
                    t_table_mapping.("Reference"){nn} = t_table_rf.PublicationRef_{kk};
                    t_table_mapping.("Comments"){nn} = t_table_rf.Comments_Notes{kk};
                    
                    
                    nn = nn + 1;
                    
                end
                
            else
                
                if strcmp(t_table3.MHCAlleleClass{mm},'-N/A-')==0
                    
                    t_table_mapping.("IEDB ID"){nn} = t_table3.IEDB(mm);
                    t_table_mapping.Protein{nn} = t_table_rf.Protein{kk};
                    t_table_mapping.("SARS Epitope"){nn} = t_table3.Epitope{mm};
%                     t_table_mapping.("SARS T Cell Assay"){nn} = t_table3.TCellAssay{mm};
                    t_table_mapping.("SARS HLA Class"){nn} = t_table3.MHCAlleleClass{mm};
                    t_table_mapping.("SARS HLA Alleles"){nn} = t_table3.MHCAlleleNames{mm};
                    t_table_mapping.("SARS Number of HLAs")(nn) = t_table3.NumMHCs(mm);
                    t_table_mapping.("SARS2 Epitope/Peptide"){nn} = t_table_rf.Peptide{kk}{:};
                    t_table_mapping.("SARS2 Host"){nn} = t_table_rf.Host{kk};
                    t_table_mapping.("SARS2 HLA Alleles"){nn} = t_table_rf.HLAAlleleInformation__Separated_{kk};
                    t_table_mapping.("SARS2 HLA Class"){nn} = t_table_rf.HLAClassInformation_I_II_{kk};
                    t_table_mapping.("SARS2 Number of Responses")(nn) = t_table_rf.No_OfSubjectsResponded(kk);
                    t_table_mapping.("SARS2 Total Subjects")(nn) = t_table_rf.No_OfSubjectsResponded(kk);
                    t_table_mapping.("RF")(nn) = t_table_rf.ResponseFreq(kk);
                    t_table_mapping.("RF CI min")(nn) = t_table_rf.ResponseFreq_CI_min(kk);
                    t_table_mapping.("RF CI max")(nn) = t_table_rf.ResponseFreq_CI_max(kk);
                    
                    t_table_mapping.("Epitope (refined)"){nn} = t_table3.Epitope{mm};
                    t_table_mapping.("doi"){nn} = t_table_rf.PublicationID_doi_{kk};
                    t_table_mapping.("Reference"){nn} = t_table_rf.PublicationRef_{kk};
                    t_table_mapping.("Comments"){nn} = t_table_rf.Comments_Notes{kk};
                    
                    
                    nn = nn + 1;
                    
                end
            end
            
        end
        
%         [kk nn]
        
    end

    if nnn==nn && ~isempty(t_table_rf.HLAAlleleInformation__Separated_{kk})
            
        t_table_mapping.("IEDB ID"){nn} = '-';
        t_table_mapping.Protein{nn} = t_table_rf.Protein{kk};
        t_table_mapping.("SARS Epitope"){nn} = '-';
        t_table_mapping.("SARS T Cell Assay"){nn} = '-';
        t_table_mapping.("SARS HLA Class"){nn} = '-';
        t_table_mapping.("SARS HLA Alleles"){nn} = '-';
        t_table_mapping.("SARS Number of HLAs")(nn) = 0;
        t_table_mapping.("SARS2 Epitope/Peptide"){nn} = t_table_rf.Peptide{kk}{:};
        t_table_mapping.("SARS2 Host"){nn} = t_table_rf.Host{kk};
        t_table_mapping.("SARS2 HLA Alleles"){nn} = t_table_rf.HLAAlleleInformation__Separated_{kk};
        t_table_mapping.("SARS2 HLA Class"){nn} = t_table_rf.HLAClassInformation_I_II_{kk};
        t_table_mapping.("SARS2 Number of Responses")(nn) = t_table_rf.No_OfSubjectsResponded(kk);
        t_table_mapping.("SARS2 Total Subjects")(nn) = t_table_rf.No_OfSubjectsResponded(kk);
        t_table_mapping.("RF")(nn) = t_table_rf.ResponseFreq(kk);
        t_table_mapping.("RF CI min")(nn) = t_table_rf.ResponseFreq_CI_min(kk);
        t_table_mapping.("RF CI max")(nn) = t_table_rf.ResponseFreq_CI_max(kk);
        
        t_table_mapping.("Epitope (refined)"){nn} = t_table_rf.Peptide{kk}{:};
        t_table_mapping.("doi"){nn} = t_table_rf.PublicationID_doi_{kk};
        t_table_mapping.("Reference"){nn} = t_table_rf.PublicationRef_{kk};
        t_table_mapping.("Comments"){nn} = t_table_rf.Comments_Notes{kk};
        
        nn = nn+1;
        
    end
    
    [kk  nn]
end


%% forming a refined list (combined from sars and sars) of HLA alleles associated with each entry in t_table_mapping

t_table_mapping.("SARS2 HLA Alleles")(find(cellfun(@isempty,t_table_mapping.("SARS2 HLA Alleles")))) = {'-'};

sars2_alleles_cell = cell(1,size(t_table_mapping,1));
sars_alleles_cell = cell(1,size(t_table_mapping,1));
alleles_refined = cell(1,size(t_table_mapping,1));
num_alleles_refined = zeros(1,size(t_table_mapping,1));

for kk = 1:size(t_table_mapping,1)
    
    sars2_alleles = t_table_mapping.("SARS2 HLA Alleles"){kk};
    sars_alleles = t_table_mapping.("SARS HLA Alleles"){kk};
    
    if strcmp(sars2_alleles,'-') == 1
        num_sars2_alleles = 0;
        sars2_alleles_cell{kk} = '-';
    else
        if sum(sars2_alleles=='/') == 0
            sars2_alleles_cell{kk}{1} = sars2_alleles;
        else
            num_sars2_alleles = sum(sars2_alleles=='/')+1;
            indices_slash = find(sars2_alleles=='/');
            sars2_alleles_cell{kk}{1} = sars2_alleles(1:indices_slash(1)-1);
            for mm = 1:num_sars2_alleles-2
                sars2_alleles_cell{kk}{mm+1} = sars2_alleles(indices_slash(mm)+1:indices_slash(mm+1)-1);
            end
            sars2_alleles_cell{kk}{end+1} = sars2_alleles(indices_slash(end)+1:end);
        end
    end
    
    if strcmp(sars_alleles,'-') == 1
        num_sars_alleles = 0;
        sars_alleles_cell{kk} = '-';
    else
        if sum(sars_alleles=='/') == 0
            sars_alleles_cell{kk}{1} = sars_alleles;
        else
            num_sars_alleles = sum(sars_alleles=='/')+1;
            indices_slash = find(sars_alleles=='/');
            sars_alleles_cell{kk}{1} = sars_alleles(1:indices_slash(1)-1);
            for mm = 1:num_sars_alleles-2
                sars_alleles_cell{kk}{mm+1} = sars_alleles(indices_slash(mm)+1:indices_slash(mm+1)-1);
            end
            sars_alleles_cell{kk}{end+1} = sars_alleles(indices_slash(end)+1:end);
        end
    end
    
    %finding unique set of alleles after combining sars and sars2 alleles
    if sars_alleles == '-'
        alleles_refined{kk} = unique(sars2_alleles_cell{kk});
    elseif sars2_alleles == '-'
        alleles_refined{kk} = unique(sars_alleles_cell{kk});
    else
        alleles_refined{kk} = unique([sars2_alleles_cell{kk} sars_alleles_cell{kk}]);
    end
    
    num_alleles_refined(kk) = length(alleles_refined{kk});
    
    %putting in t_table_mapping
    
    for nn = 1:length(alleles_refined{kk})
        if nn == 1
            t_table_mapping.("HLA Alleles (refined)"){kk} = alleles_refined{kk}{nn};   
        else
            t_table_mapping.("HLA Alleles (refined)"){kk} = ...
                sprintf('%s/%s',t_table_mapping.("HLA Alleles (refined)"){kk},alleles_refined{kk}{nn});
        end
    end
    
    t_table_mapping.("Number of HLAs (refined)")(kk) = num_alleles_refined(kk);
    
end

writetable(t_table_mapping,'Processed Data/SARS-CoV-2-T-cell-epitopes-Final-Mapped-SARS-WithHLAInfo_v2.xlsx')

%% Separating peptides excel file and extra alleles exact epitopes

indices_sars2_no_hla_info = find(strcmp(t_table_mapping.("SARS2 HLA Alleles"),'-'));
%find(cellfun(@isempty,(t_table_mapping.("SARS2 HLA Alleles"))))
t_table_immunogenic_peptides_mapped = t_table_mapping(indices_sars2_no_hla_info,:);

t_table_immunogenic_peptides_mapped.Properties.VariableNames{8} = 'SARS2 Peptide';

unique_sars_epitopes = unique(t_table_immunogenic_peptides_mapped.("SARS Epitope"));
unique_sars2_peptides = unique(t_table_immunogenic_peptides_mapped.("SARS2 Peptide"));

no_unique_sars_epitopes = length(unique_sars_epitopes)
no_unique_sars2_peptides = length(unique_sars2_peptides)

writetable(t_table_immunogenic_peptides_mapped,'Processed Data/SARS-CoV-2-T-cell-peptides-Mapped-SARS.xlsx')


%%
indices_sars2_exact_epitopes = setdiff(1:size(t_table_mapping,1),indices_sars2_no_hla_info);
t_table_exact_epitopes_mapped = t_table_mapping(indices_sars2_exact_epitopes,:);

t_table_exact_epitopes_mapped.Properties.VariableNames{8} = 'SARS2 Epitope';

indices_same_sars_sars2 = find(strcmp(t_table_exact_epitopes_mapped.("SARS Epitope"),t_table_exact_epitopes_mapped.("SARS2 Epitope")));
indices_dashes = find(strcmp(t_table_exact_epitopes_mapped.("SARS Epitope"),'-'))
% indices_notDashes = setdiff(1:size(t_table_exact_epitopes_mapped,1),indices_dashes)
% t_table_exact_epitopes_mapped2 = t_table_exact_epitopes_mapped([indices_dashes; indices_same_sars_sars2],:)


writetable(t_table_exact_epitopes_mapped,'Processed Data/SARS-CoV-2-T-cell-epitopes-Mapped-SARS.xlsx')


