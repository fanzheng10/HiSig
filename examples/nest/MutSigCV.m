%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%   MutSigCV                                                                          %
%   v1.41                                                                             %
%                                                                                     %
%   (C) 2008-2016 Mike Lawrence & Gaddy Getz                                          %
%   Broad Institute of MIT and Harvard                                                %
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% By downloading the PROGRAM you agree to the following terms of use:
%% 
%% BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
%% FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
%% 
%% This Agreement is made between the Broad Institute, Inc. with a principal 
%% address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the 
%% LICENSEE and is effective at the date the downloading is completed 
%% ("EFFECTIVE DATE").
%% WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, 
%% and BROAD wishes to have this PROGRAM utilized in the public interest, 
%% subject only to the royalty-free, nonexclusive, nontransferable license 
%% rights of the United States Government pursuant to 48 CFR 52.227-14; and
%% WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant 
%% a license on the following terms and conditions.
%% NOW, THEREFORE, in consideration of the promises and covenants made herein, 
%% the parties hereto agree as follows:
%% 
%% 1. DEFINITIONS
%% 1.1	"PROGRAM" shall mean copyright in the object code and source code 
%% known as MutSig and related documentation, if any, as they exist on the 
%% EFFECTIVE DATE and can be downloaded from 
%% http://www.broadinstitute.org/cancer/cga/MutSig on the EFFECTIVE DATE.
%% 
%% 2. LICENSE
%% 2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to 
%% LICENSEE, solely for academic non-commercial research purposes, a non-
%% exclusive, non-transferable license to: (a) download, execute and display 
%% the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
%% 
%% LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-
%% free, irrevocable license to any LICENSEE bug fixes or modifications to the 
%% PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE 
%% agrees to provide any such modifications and bug fixes to BROAD promptly 
%% upon their creation.
%% 
%% The LICENSEE may apply the PROGRAM in a pipeline to data owned by users 
%% other than the LICENSEE and provide these users the results of the PROGRAM 
%% provided LICENSEE does so for academic non-commercial purposes only.  For 
%% clarification purposes, academic sponsored research is not a commercial use 
%% under the terms of this Agreement.
%% 
%% 2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or 
%% distribute the PROGRAM, in whole or in part, without prior written 
%% permission from BROAD.  LICENSEE shall ensure that all of its users agree 
%% to the terms of this Agreement.  LICENSEE further agrees that it shall not 
%% put the PROGRAM on a network, server, or other similar technology that may 
%% be accessed by anyone other than the LICENSEE and its employees and users 
%% who have agreed to the terms of this agreement.
%% 
%% 2.3  License Limitations. Nothing in this Agreement shall be construed to 
%% confer any rights upon LICENSEE by implication, estoppel, or otherwise to 
%% any computer software, trademark, intellectual property, or patent rights 
%% of BROAD, or of any other entity, except as expressly granted herein. 
%% LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for 
%% any commercial purpose, including without limitation, as the basis of a 
%% commercial software or hardware product or to provide services. LICENSEE 
%% further agrees that the PROGRAM shall not be copied or otherwise adapted in 
%% order to circumvent the need for obtaining a license for use of the 
%% PROGRAM.  
%% 
%% 3. OWNERSHIP OF INTELLECTUAL PROPERTY 
%% LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. 
%% The PROGRAM is marked with the following BROAD copyright notice and notice 
%% of attribution to contributors. LICENSEE shall retain such notice on all 
%% copies.  LICENSEE agrees to include appropriate attribution if any results 
%% obtained from use of the PROGRAM are included in any publication.
%% 
%% Copyright 2012 Broad Institute, Inc.
%% Notice of attribution:  The MutSig program was made available through the 
%% generosity of the Cancer Genome Analysis group at the Broad Institute, Inc. 
%% 
%% LICENSEE shall not use any trademark or trade name of BROAD, or any 
%% variation, adaptation, or abbreviation, of such marks or trade names, or 
%% any names of officers, faculty, students, employees, or agents of BROAD 
%% except as states above for attribution purposes.
%% 
%% 4. INDEMNIFICATION
%% LICENSEE shall indemnify, defend, and hold harmless BROAD, and their 
%% respective officers, faculty, students, employees, associated investigators 
%% and agents, and their respective successors, heirs and assigns, 
%% ("Indemnitees"), against any liability, damage, loss, or expense (including 
%% reasonable attorneys fees and expenses) incurred by or imposed upon any of 
%% the Indemnitees in connection with any claims, suits, actions, demands or 
%% judgments arising out of any theory of liability (including, without 
%% limitation, actions in the form of tort, warranty, or strict liability and 
%% regardless of whether such action has any factual basis) pursuant to any 
%% right or license granted under this Agreement.
%% 
%% 5. NO REPRESENTATIONS OR WARRANTIES
%% THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR 
%% WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR 
%% IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, 
%% FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT 
%% OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES 
%% OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER 
%% LITERATURE MAY BE ISSUED FROM TIME TO TIME.
%% IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, 
%% AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR 
%% CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC 
%% DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD 
%% SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF 
%% THE POSSIBILITY OF THE FOREGOING.
%% 
%% 6. ASSIGNMENT
%% This Agreement is personal to LICENSEE and any rights or obligations 
%% assigned by LICENSEE without the prior written consent of BROAD shall be 
%% null and void.
%% 
%% 7. MISCELLANEOUS
%% 7.1 Export Control. LICENSEE gives assurance that it will comply with all 
%% United States export control laws and regulations controlling the export of 
%% the PROGRAM, including, without limitation, all Export Administration 
%% Regulations of the United States Department of Commerce. Among other 
%% things, these laws and regulations prohibit, or require a license for, the 
%% export of certain types of software to specified countries.
%% 7.2 Termination. LICENSEE shall have the right to terminate this Agreement 
%% for any reason upon prior written notice to BROAD. If LICENSEE breaches any 
%% provision hereunder, and fails to cure such breach within thirty (30) days, 
%% BROAD may terminate this Agreement immediately. Upon termination, LICENSEE 
%% shall provide BROAD with written assurance that the original and all copies 
%% of the PROGRAM have been destroyed, except that, upon prior written 
%% authorization from BROAD, LICENSEE may retain a copy for archive purposes.
%% 7.3 Survival. The following provisions shall survive the expiration or 
%% termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 
%% 7.3, and 7.4.
%% 7.4 Notice. Any notices under this Agreement shall be in writing, shall 
%% specifically refer to this Agreement, and shall be sent by hand, recognized 
%% national overnight courier, confirmed facsimile transmission, confirmed 
%% electronic mail, or registered or certified mail, postage prepaid, return 
%% receipt requested.  All notices under this Agreement shall be deemed 
%% effective upon receipt. 
%% 7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, 
%% supplemented, or otherwise modified only by means of a written instrument 
%% signed by all parties. Any waiver of any rights or failure to act in a 
%% specific instance shall relate only to such instance and shall not be 
%% construed as an agreement to waive any rights or fail to act in any other 
%% instance, whether or not similar. This Agreement constitutes the entire 
%% agreement among the parties with respect to its subject matter and 
%% supersedes prior agreements or understandings between the parties relating 
%% to its subject matter. 
%% 7.6 Binding Effect; Headings. This Agreement shall be binding upon and 
%% inure to the benefit of the parties and their respective permitted 
%% successors and assigns. All headings are for convenience only and shall not 
%% affect the meaning of any provision of this Agreement.
%% 7.7 Governing Law. This Agreement shall be construed, governed, interpreted 
%% and applied in accordance with the internal laws of the Commonwealth of 
%% Massachusetts, U.S.A., without regard to conflict of laws principles.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MutSigCV(mutation_file,coverage_file,covariate_file,output_filestem,varargin)

  params1 = 'mutation_file,coverage_file,covariate_file,output_filestem';
  params2 = 'mutation_type_dictionary_file,chr_files_directory,categ_flag';
  if nargin<4, error(['usage: MutSig_preprocess(' params1 ')']); end
  if nargin>7, error(['usage: MutSig_preprocess(' params1 ',' params2 ')']); end

  fprintf('\n');
  fprintf('======================================\n');
  fprintf('  MutSigCV\n');
  fprintf('  v1.4\n\n');
  fprintf('  (c) Mike Lawrence and Gaddy Getz\n');
  fprintf('  Broad Institute of MIT and Harvard\n');
  fprintf('======================================\n');

  % (1) PREPROCESS
  fprintf('\n\n');
  fprintf('MutSigCV: PREPROCESS\n');
  fprintf('--------------------\n');

  MutSig_preprocess(mutation_file,coverage_file,covariate_file,output_filestem,varargin{:})

  % (2) RUN
  fprintf('\n\n');
  fprintf('MutSigCV: RUN\n');
  fprintf('-------------\n');

  mutation_file = [output_filestem '.mutations.txt'];
  coverage_file = [output_filestem '.coverage.txt'];
  output_file =   [output_filestem '.sig_genes.txt'];
  MutSig_runCV(mutation_file,coverage_file,covariate_file,output_file);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MutSig_preprocess(mutation_file,coverage_file,covariate_file,output_filestem,mutation_type_dictionary_file,chr_files_directory,categ_flag)

  params1 = 'mutation_file,coverage_file,covariate_file,output_filestem';
  params2 = 'mutation_type_dictionary_file,chr_files_directory,categ_flag';
  if nargin<4, error(['usage: MutSig_preprocess(' params1 ')']); end
  if nargin>7, error(['usage: MutSig_preprocess(' params1 ',' params2 ')']); end

  if ~exist('mutation_type_dictionary_file','var'), mutation_type_dictionary_file = 'mutation_type_dictionary_file.txt'; end
  if ~exist('chr_files_directory','var'), chr_files_directory = 'chr_files'; end
  if ~exist('categ_flag','var'), categ_flag = nan; end
  if ischar(categ_flag)
    categ_flag = str2double(categ_flag);
    if isnan(categ_flag), error('categ_flag should be numeric'); end
  end

  % first, ensure output_filestem directory exists and can be written to
  [outpath] = fileparts(output_filestem);
  if ~isempty(outpath) && ~exist(outpath,'dir'), mkdir(outpath); end
  ensure_writeable([output_filestem '.test.txt']);

  % load MUTATION FILE and make sure has gene+patient
  
  fprintf('Loading mutation_file...\n');
  M = load_struct(mutation_file);

  % GENE
  if isfield(M,'gene') && isfield(M,'Hugo_Symbol')
    fprintf('NOTE:  Both "gene" and "Hugo_Symbol" are present in mutation_file.  Using "gene".\n');
  elseif isfield(M,'gene')
    % OK
  elseif isfield(M,'Hugo_Symbol')
    M.gene = M.Hugo_Symbol;
  else
    error('mutation_file lacks "gene" or "Hugo_Symbol" column.');
  end
  
  % PATIENT
  if isfield(M,'patient') && isfield(M,'Tumor_Sample_Barcode')
    fprintf('NOTE:  Both "patient" and "Tumor_Sample_Barcode" are present in mutation_file.  Using "patient".\n');
  elseif isfield(M,'patient')
    % OK
  elseif isfield(M,'Tumor_Sample_Barcode')
    M.patient = M.Tumor_Sample_Barcode;
  else
    error('mutation_file lacks "patient" or "Tumor_Sample_Barcode" column.');
  end
  if length(unique(M.patient))<2, error('MutSig is not applicable to single patients.\n'); end

  % LOAD COVERAGE FILE and make sure has gene+effect+categ

  fprintf('Loading coverage file...\n');
  C = load_struct_specify_string_cols(coverage_file,1:3); % gene effect categ are all strings
  if ~isfield(C,'gene'), error('no "gene" column in coverage_file'); end
  if ~isfield(C,'effect') && isfield(C,'zone'), C = rename_field(C,'zone','effect'); end
  if ~isfield(C,'effect'), error('no "effect" column in coverage_file'); end
  C.effect = regexprep(C.effect,'^flank.*','noncoding');
  if any(~ismember(unique(C.effect),{'noncoding','silent','nonsilent'}))
    error('in coverage_file, "effect" must be one of noncoding/silent/nonsilent');
  end
  if ~isfield(C,'categ'), error('no "categ" column in coverage_file'); end
  f = fieldnames(C); coverage_patient_names = f(4:end);
  
  % MAKE SURE COVARIATES FILE EXISTS
  if ~exist(covariate_file,'file'), error('covariates_file not found'); end

  %%%%%%%%%%%%
  % EFFECT   %
  %%%%%%%%%%%%
  
  fprintf('Processing mutation "effect"...\n');

  if isfield(M,'is_coding') || isfield(M,'is_silent')
    fprintf('NOTE:  This version now ignores "is_coding" and "is_silent".\n');
    fprintf('       Requires Variant_Classification/type column and mutation_type_dictionary so we can assign nulls.\n');
    M = rmfield_if_exist(M,{'is_coding','is_silent'});
  end
  if isfield(M,'effect')
    fprintf('Will use the pre-existing "effect" column.\n');
    M.effect = regexprep(M.effect,'^flank.*','noncoding');
    if any(~ismember(unique(M.effect),{'noncoding','silent','nonsilent','null'}))
      error('in mutation_file, "effect" must be one of noncoding/silent/nonsilent/null.');
    end
  else
    % try to calculate EFFECT column from Variant Classification/Type
    if ~isfield(M,'Variant_Classification') && isfield(M,'type')
      M.Variant_Classification = M.type;
    end
    if ~isfield(M,'Variant_Classification')
      error('mutation_file is missing Variant_Classification');
    end
    if isempty(mutation_type_dictionary_file) || ~exist(mutation_type_dictionary_file,'file')
      error('missing mutation_type_dictionary_file');
    end
    dict = load_struct(mutation_type_dictionary_file);
    require_fields(dict,{'Variant_Classification','effect'});
    M.effect = mapacross(upper(M.Variant_Classification),upper(dict.Variant_Classification),dict.effect,'unknown');
    bad = find(strcmp('unknown',M.effect));
    if length(bad)>0
      fprintf('WARNING:  %d/%d mutations could not be mapped to effect using mutation_type_dictionary_file:\n',length(bad),slength(M));
      count(M.Variant_Classification(bad),1);
      fprintf('          They will be removed from the analysis.\n');
      M = reorder_struct_exclude(M,bad);
    end
    if slength(M)==0, error('No mutations left!'); end
  end

  %%%%%%%%%%%%
  % CATEG    %
  %%%%%%%%%%%%
  
  fprintf('Processing mutation "categ"...\n');
  
  % see if coverage file is on the FULL192 set
  ucc = unique(C.categ);
  names192 = generate_192_categ_names();
  coverage_is_on_full192 = (length(ucc)==192 && length(intersect(ucc,names192))==192); 
  categs_already_present = false;
  
  if isfield(M,'categ')
    % categ already provided in M:
    % make sure categories in C are same as in M
    mcc = unique(M.categ);
    ucc = unique(C.categ);
    if length(ucc)~=length(mcc) || any(~ismember(ucc,mcc)) || any(~ismember(mcc,ucc))
      fprintf('"categ" of mutation_file does not match coverage_file.  Ignoring it.\n');
      M = rmfield(M,'categ');
    else
      categs_already_present = true;
    end
  end

  % find out if we can do category discovery
  can_do_category_discovery = true;
  if ~coverage_is_on_full192, can_do_category_discovery = false; reason = 'coverage_file not on full192'; end
  if ~isfield(M,'chr') && isfield(M,'Chromosome'), M.chr=M.Chromosome; end
  if ~isfield(M,'start') && isfield(M,'Start_position'), M.start = M.Start_position; end
  if ~isfield(M,'start') && isfield(M,'Start_Position'), M.start = M.Start_Position; end
  if ~isfield(M,'ref_allele') && isfield(M,'Reference_Allele'), M.ref_allele = M.Reference_Allele; end
  if ~isfield(M,'newbase')
    if isfield(M,'Tumor_Seq_Allele1')
      M.newbase = M.Tumor_Seq_Allele1;
      if isfield(M,'Tumor_Seq_Allele2')
        idx = find(strcmp(M.ref_allele,M.Tumor_Seq_Allele1));
        M.newbase(idx) = M.Tumor_Seq_Allele2(idx);
      end
    end
  end
  if ~isfield(M,'chr') || ~isfield(M,'start') || ~isfield(M,'ref_allele') || ~isfield(M,'newbase')
    can_do_category_discovery = false;
    reason = 'Chromosome/Start_position/Reference_Allele/Tumor_Seq_Allele1/Tumor_Seq_Allele2 missing from mutation_file';
  else
    % make "start" numeric
    M = make_numeric(M,'start');
    bad = find(isnan(M.start));
    if length(bad)>0
      fprintf('WARNING:  %d/%d mutations had non-numeric Start_position.  Excluding them from analysis.\n',length(bad),slength(M));
      M = reorder_struct_exclude(M,bad);
    end
    if slength(M)==0, error('No mutations left!\n'); end
    % see if chromosome files are available
    if isempty(chr_files_directory)
      can_do_category_discovery = false;
      reason = 'no chr_files_directory available';
    else
      f1 = direc([chr_files_directory '/chr*.txt']);
      [uchr tmp M.chr_idx] = unique(M.chr);
      uchr = regexprep(uchr,'^chr','');
      f2 = regexprep(uchr,'^(.*)$',[chr_files_directory '/chr$1.txt']);
      chr_file_available = ismember(f2,f1);
      if ~any(chr_file_available)
        can_do_category_discovery = false;
        reason = 'no chr_files available';
      else
        % remove mutations on weird chromosomes
        bad = find(~chr_file_available(M.chr_idx));
        if length(bad)>0
          fprintf('WARNING:  %d/%d mutations are on chromosomes not found in chr_files_directory.  Excluding them from analysis.\n',length(bad),slength(M));
          M = reorder_struct_exclude(M,bad);
        end
        if slength(M)==0, error('No mutations left!\n'); end
      end
    end
  end

  % DECIDE WHAT TO DO ABOUT CATEGORIES
  % METHODS:   1. use the existing categories
  %            2. have only one category for non-nulls
  %            3. discover k categories for non-nulls
  
  if categ_flag==0
    % requested use of categories already present
    if ~categs_already_present
      error('when setting categ_flag==0, "categ" column must be already present in mutation_file.');
    end
    method = 1;
  elseif categ_flag==1
    % requested only one category for non-nulls
    method = 2;
  elseif categ_flag>1
    % requested category discovery
    if categ_flag>6, fprintf('NOTE:  maximum categories that can be discovered is 6.\n'); categ_flag=6; end
    if ~can_do_category_discovery
      error('unable to perform category discovery, because %s.',reason);
    end
  elseif isnan(categ_flag)
    % no method specified: let's choose what to do
    if can_do_category_discovery
      % let's do it.  discover four categories by default.
      method = 3;
      ncategs = 4;
    else
      if categs_already_present
        method = 1; % use them
      else
        fprintf('NOTE:  unable to perform category discovery, because %s.',reason);
        method = 2; % will put all non-null in one category "missense"
      end
    end
  end
  
  if method==1
    fprintf('Will use the categories already present.\n');
    K=[];
    K.name = unique(M.categ);
    
  elseif method==2
    fprintf('Will use two categories: missense and null+indel.\n');
    
    K1 = [];
    K1.left = {'ACGT'};
    K1.from = {'AC'};
    K1.change = {'in'};
    K1.right = {'ACGT'};
    K1.autoname = {'missense'};
    K1.name = {'missense'};
    K1.type = {'point'};
    
    K2 = [];
    K2.left = {'ACGT'};
    K2.from = {'AC'};
    K2.change = {'in'};
    K2.right = {'ACGT'};
    K2.autoname = {'null+indel'};
    K2.name = {'null+indel'};
    K2.type = {'non-point'};
    
    K = concat_structs_keep_all_fields({K1,K2});
    
    % assign categories
    M.categ = nansub(K.name,1 + strcmp('null',M.effect));
    
    % collapse coverage
    fprintf('Collapsing coverage...\n');
    [tmp tmp C.categ_idx] = unique(C.categ);
    C = sort_struct(C,{'gene','effect','categ_idx'});
    ug = unique(C.gene); ng = length(ug);
    ue = unique(C.effect); ne = length(ue);
    nk = slength(K);
    
    idx = find(C.categ_idx<=nk);
    C2 = reorder_struct(C,idx);
    C2.categ = nansub(K.name,C2.categ_idx);
    C2 = keep_fields(C2,{'gene','effect','categ'});
    
    np = length(coverage_patient_names);
    for p=1:np
      oldcov = reshape(C.(coverage_patient_names{p}),[192 ne ng]);
      newcov = repmat(sum(oldcov,1),2,1);               % both = total territory
      C2.(coverage_patient_names{p}) = newcov(:);
    end
    C=C2; clear C2;
    
  elseif method==3
    
    %%%%%%%%%%%%%%%%%%%%%%
    % CATEGORY DISCOVERY %
    %%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('Preparing for category discovery...\n');
    
    % STEP1
    % COVERAGE: get total on context65
    
    % get total coverage
    C.totcov = C.(coverage_patient_names{1});
    for i=2:length(coverage_patient_names)
      C.totcov = C.totcov + C.(coverage_patient_names{i});
    end
    npm = length(unique(M.patient));
    if length(coverage_patient_names)==1 && npm>1
      C.totcov = C.totcov * npm;
    end
    
    % will use only the coding mutations+coverage to do this
    C.is_coding = (strcmp(C.effect,'nonsilent')|strcmp(C.effect,'silent'));
    
    % collapse coverage to 192
    X = [];
    [X.categ tmp C.categ_idx] = unique(C.categ);
    X = parse_in(X,'categ','(.).(.)..(.).(.)',{'left','from','to','right'});
    X.yname = cell(slength(X),1);
    X.N = nan(slength(X),1);
    for i=1:slength(X)
      X.yname{i} = [X.from{i} ' in ' X.left{i} '_' X.right{i}];
      X.N(i) = sum(C.totcov(C.is_coding & C.categ_idx==i));
    end
    X.newbase_idx = listmap(X.to,{'A','C','G','T'});
    
    % collapse coverage to context65
    Y = generate_categ_context65_names();
    X.context65 = listmap(X.yname,Y.name);
    Y.N = nan(slength(Y),1);
    for i=1:slength(Y)
      Y.N(i) = sum(X.N(X.context65==i));
    end
    Y.N = round(Y.N/3);
    
    % STEP2
    % MUTATIONS: get context65 by looking up from reference genome
    fprintf('Looking up trinucleotide contexts from chr_files...\n');
    M.triplet = repmat({'---'},slength(M),1);
    for ci=1:length(uchr), fprintf(' %d/%d',ci,length(uchr));
      midx = find(M.chr_idx==ci);
      chrfile = f2{ci};
      d = dir(chrfile);
      if isempty(d), continue; end
      filesize = d.bytes;
      ff = fopen(chrfile);
      for i=1:length(midx),mi=midx(i);
        leftpos = M.start(mi)-1;
        if leftpos>1 && leftpos+2<=filesize
          status = fseek(ff, leftpos-1, 'bof');
          M.triplet{mi} = fgets(ff,3);
        end
      end
    end, fprintf('\n');
    M.triplet = upper(M.triplet);
    M = parse_in(M,'triplet','^.(.).$','triplet_middle');
    midx = find(~strcmp('-',M.ref_allele)&~strcmp('-',M.newbase));
    matchfrac = mean(strcmpi(M.ref_allele(midx),M.triplet_middle(midx)));
    if matchfrac<0.9
      adj = 'possible'; if matchfrac<0.7, adj = 'probable'; end
      error('%s build mismatch between mutation_file and chr_files', adj);
    end
    M.yname = regexprep(M.triplet,'^(.)(.)(.)$','$2 in $1_$3');
    M.context65 = listmap(M.yname,Y.name);
    M.newbase_idx = listmap(regexprep(M.newbase,'^(.).*','$1'),{'A','C','G','T'});
    midx = find(~strcmp('-',M.ref_allele) & ~strcmp('-',M.newbase) & M.context65>=1 & M.context65<=65 & M.newbase_idx>=1 & M.newbase_idx<=4);
    Y.n = hist2d(M.context65(midx),M.newbase_idx(midx),1:65,1:4);
    
    % STEP3
    % Category Discovery
    Nn = collapse_Nn_65_to_32([Y.N Y.n]);
    PP=[]; PP.max_k = ncategs;
    PP.mutcategs_report_filename = [output_filestem '.mutcateg_discovery.txt'];
    Ks = find_mut_categs(Nn,PP);
    K = Ks{ncategs};
    c = assign_65x4_to_categ_set(K);
    X.kidx = nan(slength(X),1);
    for i=1:slength(X)
      X.kidx(i) = find(c(X.context65(i),:,X.newbase_idx(i)));
    end
    
    % STEP4 
    % assign mutation categories
    fprintf('Assigning mutation categories...\n');
    M.categ = repmat({'---'},slength(M),1);
    for i=1:slength(X)
      idx=find(M.context65==X.context65(i) & M.newbase_idx==X.newbase_idx(i));
      M.categ(idx) = repmat(K.name(X.kidx(i)),length(idx),1);
    end
    % add null+indel category
    K2 = [];
    K2.left = {'ACGT'};
    K2.from = {'AC'};
    K2.change = {'in'};
    K2.right = {'ACGT'};
    K2.autoname = {'null+indel'};
    K2.name = {'null+indel'};
    K2.type = {'non-point'};
    K = concat_structs_keep_all_fields({K,K2});
    K.N(end) = sum(K.N(1:end-1));
    midx = find(strcmp('null',M.effect));
    M.categ(midx) = repmat(K2.name(end),length(midx),1);
    K.n(end) = length(midx);
    K.rate(end) = K.n(end)/K.N(end);
    K.relrate(end) = K.rate(end)/K.rate(1)*K.relrate(1);
    
    % STEP5
    % collapse coverage
    
    fprintf('Collapsing coverage...\n');
    C = sort_struct(C,{'gene','effect','categ_idx'});
    ug = unique(C.gene); ng = length(ug);
    ue = unique(C.effect); ne = length(ue);
    nk = slength(K);
    
    idx = find(C.categ_idx<=nk);
    C2 = reorder_struct(C,idx);
    C2.categ = nansub(K.name,C2.categ_idx);
    C2 = keep_fields(C2,{'gene','effect','categ'});
    
    np = length(coverage_patient_names);
    for p=1:np
      oldcov = reshape(C.(coverage_patient_names{p}),[192 ne ng]);
      newcov = nan(nk,ne,ng);
      for ki=1:nk
        if ki==nk
          cidx = 1:192;  % null+indel = total territory
        else
          cidx = find(X.kidx==ki);
        end
        newcov(ki,:,:) = sum(oldcov(cidx,:,:),1);
      end
      C2.(coverage_patient_names{p}) = newcov(:);
    end
    C=C2; clear C2;
    
  end
  
  % SAVE OUTPUT FILES
  
  fprintf('Writing preprocessed files.\n');
  
  % (1) mutation file
  M = rmfield_if_exist(M,{'newbase','chr_idx','triplet','yname','context65','newbase_idx','context65_right','triplet_middle'});
  save_struct(M,[output_filestem '.mutations.txt']);
  
  % (2) coverage file
  save_struct(C,[output_filestem '.coverage.txt']);
  
  % (3) categories file
  save_struct(K,[output_filestem '.categs.txt']);
  
  fprintf('MutSig_preprocess finished.\n');
  
end % end of MutSig_preprocess

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

function MutSig_runCV(mutation_file,coverage_file,covariate_file,output_file)

  if nargin~=4, error('usage: MutSig_runCV(mutation_file,coverage_file,covariate_file,output_file)'); end

  % first, ensure output_file directory exists
  [outpath] = fileparts(output_file);
  if ~isempty(outpath) && ~exist(outpath,'dir'), mkdir(outpath); end
  
  % load MUTATION FILE and make sure all required fields are present

  fprintf('Loading mutation_file...\n');
  M = load_struct(mutation_file);
  
  % GENE
  if isfield(M,'gene') && isfield(M,'Hugo_Symbol')
    fprintf('NOTE:  Both "gene" and "Hugo_Symbol" are present in mutation_file.  Using "gene".\n');
  elseif isfield(M,'gene')
    % OK
  elseif isfield(M,'Hugo_Symbol')
    M.gene = M.Hugo_Symbol;
  else
    error('mutation_file lacks "gene" or "Hugo_Symbol" column.');
  end
  
  % PATIENT
  if isfield(M,'patient') && isfield(M,'Tumor_Sample_Barcode')
    fprintf('NOTE:  Both "patient" and "Tumor_Sample_Barcode" are present in mutation_file.  Using "patient".\n');
  elseif isfield(M,'patient')
    % OK
  elseif isfield(M,'Tumor_Sample_Barcode')
    M.patient = M.Tumor_Sample_Barcode;
  else
    error('mutation_file lacks "patient" or "Tumor_Sample_Barcode" column.');
  end
  
  % EFFECT
  if isfield(M,'is_coding') || isfield(M,'is_silent')
    fprintf('NOTE:  This version now ignores "is_coding" and "is_silent".  Requires "effect" column.\n');
    M = rmfield_if_exist(M,{'is_coding','is_silent'});
  end
  if isfield(M,'effect')
    M.effect = regexprep(M.effect,'^flank.*','noncoding');
    if any(~ismember(unique(M.effect),{'noncoding','silent','nonsilent','null'}))
      error('in mutation_file, "effect" must be one of noncoding/silent/nonsilent/null');
    end
  else
    error('mutation_file lacks "effect" column.');
  end
  
  % CATEG
  if ~isfield(M,'categ'), error('mutation_file lacks "categ" column.'); end
  
  % LOAD COVERAGE FILE and COVARIATE FILE
  fprintf('Loading coverage file...\n');
  C = load_struct_specify_string_cols(coverage_file,1:3); % gene effect categ are all strings
  G=[]; G.gene = unique(C.gene); ng=slength(G);
  fprintf('Loading covariate file...\n');
  V = load_struct_specify_string_cols(covariate_file,1);  % gene is string
  f = fieldnames(V); cvnames = f(2:end); nv = length(cvnames);
  gidx = listmap(G.gene,V.gene);
  for i=1:length(f)
    if strcmp(f{i},'gene'), continue; end
    G.(f{i}) = nansub(V.(f{i}),gidx);
  end
  
  % make sure coverage file has required fields
  if ~isfield(C,'gene'), error('no "gene" column in coverage_file'); end
  if ~isfield(C,'effect') && isfield(C,'zone'), C = rename_field(C,'zone','effect'); end
  if ~isfield(C,'effect'), error('no "effect" column in coverage_file'); end
  C.effect = regexprep(C.effect,'^flank.*','noncoding');
  if any(~ismember(unique(C.effect),{'noncoding','silent','nonsilent'}))
    error('in coverage_file, "effect" must be one of noncoding/silent/nonsilent');
  end
  if ~isfield(C,'categ'), error('no "categ" column in coverage_file'); end
  f = fieldnames(C); coverage_patient_names = f(4:end);
  
  % remove any genes that we don't have coverage for
  badgene = setdiff(M.gene,C.gene);
  if ~isempty(badgene)
    fprintf('NOTE:  %d/%d gene names could not be mapped to coverage information.  Excluding them.\n',length(badgene),length(unique(M.gene)));
    M = reorder_struct_exclude(M,ismember(M.gene,badgene));
  end

  % make sure categories in C are same as in M
  bad = find(~ismember(M.categ,C.categ));
  if ~isempty(bad)
    fprintf('NOTE:  %d/%d mutations were outside the category set.  Excluding them.\n',length(bad),slength(M));
    if isfield(M,'ref_allele') && isfield(M,'newbase') && isfield(M,'type')
      is_probably_indel = strcmp('-',M.ref_allele(bad)) | strcmp('-',M.newbase(bad));
      is_probably_noncoding = grepmi('intron|utr|igr|flank',M.type(bad));
      nbad2 = bad(is_probably_indel&is_probably_noncoding);
      if nbad>0
        fprintf('           (%d of them are noncoding indels.)\n',length(bad2));
      end
    end
    M = reorder_struct_exclude(M,bad);
  end
  if slength(M)==0, error('No mutations left!\n'); end

  % map categories
  [K.name tmp C.categ_idx] = unique(C.categ);
  M.categ_idx = listmap(M.categ,K.name);
  ncat = slength(K);

  % make sure there is a null+indel category
  knum = str2double(K.name);
  if all(~isnan(knum))
    % all category names are numbers: assume null+indel is the last one
    % (this maintains compatibility with the LUSC test data files)
    [tmp null_categ] = max(knum);
  else
    null_categ = grepi('null|indel',K.name,1);
    if length(null_categ)==0, error('ERROR: no null/indel category.\n'); end
    if length(null_categ)>1, error('ERROR: multiple null/indel categories.\n'); end
  end

  % make sure C is sorted by the same gene order as in G
  C.gene_idx = listmap(C.gene,G.gene);
  C = sort_struct(C,'gene_idx');
  
  % map genes
  M.gene_idx = listmap(M.gene,G.gene);
  
  % regularize the sample name in the mutation table
  namebefore = M.patient;
  M.patient = regexprep(M.patient,'-Tumor$','');
  if any(~strcmp(namebefore,M.patient)), fprintf('NOTE:  Trimming "-Tumor" from patient names.\n'); end
  namebefore = M.patient;
  M.patient = regexprep(M.patient,'-','_');
  if any(~strcmp(namebefore,M.patient)), fprintf('NOTE:  Converting "-" to "_" in patient names.\n'); end

  pat=[]; [pat.name tmp M.patient_idx] = unique(M.patient);
  pat.cov_idx = listmap(pat.name,coverage_patient_names);
  
  % Fan
  fid = fopen(strrep(output_file, 'sig_genes.txt', 'patients.txt'), 'w');
  fprintf(fid, '%s\n', pat.name{:});
  fclose(fid);
  
  np = slength(pat);
  if np<2, error('MutSig is not applicable to single patients.\n'); end

  % is generic coverage data given?
  generic_column_name = 'coverage';
  if length(coverage_patient_names)>1 || (length(coverage_patient_names)==1 && ~strcmpi(coverage_patient_names{1},generic_column_name))
    if any(strcmp(coverage_patient_names,generic_column_name)), error('reserved name "%s" cannot appear in list of patient names',generic_column_name); end
    if length(coverage_patient_names)~=length(unique(coverage_patient_names)), error('patient names in coverage_file must be unique'); end
    % make sure all patients are accounted for in coverage file
    if any(isnan(pat.cov_idx)), error('some patients in mutation_file are not accounted for in coverage_file'); end
    generic_coverage_flag = false;
  else
    % we're using generic coverage
    pat.cov_idx(:) = 1;
    generic_coverage_flag = true;
  end
  
  % BUILD n and N tables
  fprintf('Building n and N tables...\n');
  
  midx = strcmpi(M.effect,'silent');
  n_silent = hist3d(M.gene_idx(midx),M.categ_idx(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
  midx = strcmpi(M.effect,'nonsilent') | strcmpi(M.effect,'null');
  n_nonsilent = hist3d(M.gene_idx(midx),M.categ_idx(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
  midx = strcmpi(M.effect,'noncoding');
  n_noncoding = hist3d(M.gene_idx(midx),M.categ_idx(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
  
  N_silent = nan(ng,ncat,np);
  N_nonsilent = nan(ng,ncat,np);
  N_noncoding = nan(ng,ncat,np);
  
  for ci=1:ncat
    silent_idx = strcmpi(C.effect,'silent') & C.categ_idx==ci;
    nonsilent_idx = strcmpi(C.effect,'nonsilent') & C.categ_idx==ci;
    noncoding_idx = strcmpi(C.effect,'noncoding') & C.categ_idx==ci;
    for pi=1:np
      cpfld = coverage_patient_names{pat.cov_idx(pi)};
      N_silent(:,ci,pi) = C.(cpfld)(silent_idx);
      N_nonsilent(:,ci,pi) = C.(cpfld)(nonsilent_idx);
      N_noncoding(:,ci,pi) = C.(cpfld)(noncoding_idx);
    end
  end

  % MAKE SURE ALL NUMBERS ARE INTEGERS
  n_silent = round(n_silent);
  n_nonsilent = round(n_nonsilent);
  n_noncoding = round(n_noncoding);
  N_silent = round(N_silent);
  N_nonsilent = round(N_nonsilent);
  N_noncoding = round(N_noncoding);

  % REMOVE MUTATIONS IN BINS WITH EXTREMELY LOW COVERAGE
  n_silent(n_silent>N_silent) = 0;
  n_nonsilent(n_nonsilent>N_nonsilent) = 0;
  n_noncoding(n_noncoding>N_noncoding) = 0;
  
  % SANITY CHECKS ON TOTALS
  tot_n_nonsilent = fullsum(n_nonsilent);
  tot_N_nonsilent = fullsum(N_nonsilent);
  tot_n_silent = fullsum(n_silent);
  tot_N_silent = fullsum(N_silent);
  tot_n_noncoding = fullsum(n_noncoding);
  tot_N_noncoding = fullsum(N_noncoding);
  tot_rate_nonsilent = tot_n_nonsilent/tot_N_nonsilent;
  tot_rate_silent = tot_n_silent/tot_N_silent;
  tot_rate_noncoding = tot_n_noncoding/tot_N_noncoding;
  tot_rate_coding = (tot_n_nonsilent+tot_n_silent)/(tot_N_nonsilent+tot_N_silent);
  
  min_tot_n_nonsilent = 50;
  min_tot_n_silent = 50;
  min_tot_n_noncoding = 50;
  min_rate_nonsilent = 1e-9;
  max_rate_nonsilent = 1e-3;
  min_rate_silent = 1e-9;
  max_rate_silent = 1e-3;
  min_rate_noncoding = 1e-9;
  max_rate_noncoding = 1e-3;
  max_abs_log2_difference_nonsilent_silent = 1.0;
  max_abs_log2_difference_noncoding_coding = 1.0;
  
  % see if silent and nonsilent are OK: if not, give warning
  if tot_n_nonsilent<min_tot_n_nonsilent || tot_n_silent<min_tot_n_silent, error('not enough mutations to analyze'); end
  if tot_rate_nonsilent<min_rate_nonsilent || tot_rate_nonsilent>max_rate_nonsilent, error('nonsilent mutation rate out of range'); end
  if tot_rate_silent<min_rate_silent || tot_rate_silent>max_rate_silent, error('silent mutation rate out of range'); end
  abs_log2_difference_nonsilent_silent = abs(log2(tot_rate_nonsilent/tot_rate_silent));
  if abs_log2_difference_nonsilent_silent>max_abs_log2_difference_nonsilent_silent, fprintf('Warning: silent and nonsilent rates are too different.\n'); end
  
  % see if noncoding is OK: if not, give warning and zero it all out
  ok = false;
  if tot_n_noncoding==0
    fprintf('NOTE:  no noncoding mutations.\n');
  else
    if tot_n_noncoding<min_tot_n_noncoding
      fprintf('WARNING:  not enough noncoding mutations to analyze\n');
    else
      if tot_rate_noncoding<min_rate_noncoding || tot_rate_noncoding>max_rate_noncoding
        fprintf('WARNING:  noncoding mutation rate out of range\n');
      else
        abs_log2_difference_noncoding_coding = abs(log2(tot_rate_noncoding/tot_rate_coding));
        if abs_log2_difference_noncoding_coding>max_abs_log2_difference_noncoding_coding
          fprintf('WARNING:  coding and noncoding rates are too different\n');
        else
          ok = true;
  end,end,end,end
  if ~ok
    fprintf('Zeroing out all noncoding mutations and coverage for the rest of the calculation.\n');
    n_noncoding(:) = 0;
    N_noncoding(:) = 0;
  end
  
  % add total columns
  n_silent(:,end+1,:) = sum(n_silent,2);
  n_nonsilent(:,end+1,:) = sum(n_nonsilent,2);
  n_noncoding(:,end+1,:) = sum(n_noncoding,2);
  N_silent(:,end+1,:) = N_silent(:,null_categ,:);          % copy total coverage from null+indel coverage
  N_nonsilent(:,end+1,:) = N_nonsilent(:,null_categ,:);
  N_noncoding(:,end+1,:) = N_noncoding(:,null_categ,:);
  
  % total across patients, save in G
  G.N_nonsilent = sum(N_nonsilent(:,end,:),3);
  G.N_silent = sum(N_silent(:,end,:),3);
  G.N_noncoding = sum(N_noncoding(:,end,:),3);
  G.n_nonsilent = sum(n_nonsilent(:,end,:),3);
  G.n_silent = sum(n_silent(:,end,:),3);
  G.n_noncoding = sum(n_noncoding(:,end,:),3);
  
  % PROCESS COVARIATES
  
  fprintf('Processing covariates...\n');
  
  V = nan(ng,nv);
  for vi=1:nv, V(:,vi) = G.(cvnames{vi}); end
  
  % convert covariate raw values to Z-scores
  Z = nan(ng,nv);
  for vi=1:nv
    missing = isnan(V(:,vi)) | isinf(V(:,vi));
    mn = mean(V(~missing,vi));
    sd = std(V(~missing,vi));
    Z(~missing,vi) = (V(~missing,vi)-mn)./sd;
  end
  
  % FIND BAGELS
  
  fprintf('Finding bagels...  ');
  
  max_neighbors = 50;
  qual_min = 0.05;
  
  G.nnei = nan(ng,1); G.x = nan(ng,1); G.X = nan(ng,1);
  
  for gi=1:ng, if ~mod(gi,1000), fprintf('%d/%d ',gi,ng); end
    
    % calculate distances from this gene
    df2 = bsxfun(@minus,Z,Z(gi,:)).^2;
    dist2 = nansum(df2,2)./sum(~isnan(df2),2);
    [tmp,ord] = sort(dist2); ord = [gi;ord(ord~=gi)];
    
    % expand bagel outward until quality falls below qual_min
    nfit=0; Nfit=0;
    for ni=0:max_neighbors, gidx = ord(ni+1);
      
      ngene = G.n_silent(gidx) + G.n_noncoding(gidx);
      Ngene = G.N_silent(gidx) + G.N_noncoding(gidx);
      if ni==0, ngene0=ngene; Ngene0=Ngene; end
      nfit=nfit+ngene; Nfit=Nfit+Ngene;
      
      % compare the gene being added to the central gene
      hc = hyge2cdf(ngene,Ngene,ngene0,Ngene0);
      qual_left = min(hc, 1-hc);
      qual = 2*qual_left;
      
      % stopping criterion: stop if this gene would drop quality below qual_min
      if ni>0 && qual<qual_min, break; end
      
      % update gene's statistics
      G.nnei(gi) = ni; G.x(gi) = nfit; G.X(gi) = Nfit;
      
    end % next neighborhood size
  end, fprintf('\n'); % next gene
  
  fprintf('Expanding to (x,X)_gcp...\n');
  
  n_gcp = n_nonsilent + n_silent + n_noncoding;
  N_gcp = N_nonsilent + N_silent + N_noncoding;
  
  n_cp = sum(n_gcp,1);
  N_cp = sum(N_gcp,1);
  
  n_c = sum(n_cp,3);
  N_c = sum(N_cp,3);
  mu_c = n_c./N_c;
  
  n_tot = n_c(end);
  N_tot = N_c(end);
  mu_tot = n_tot/N_tot;
  f_c = mu_c/mu_tot;
  f_Nc = N_c/N_tot;
  
  n_p = n_cp(:,end,:);
  N_p = N_cp(:,end,:);
  mu_p = n_p./N_p;
  f_p = mu_p/mu_tot;
  f_Np = N_p/mean(N_p);
  
  x_gcp = repmat(G.x,[1 ncat+1 np]); X_gcp = repmat(G.X,[1 ncat+1 np]);       % last column = total
  x_gcp = bsxfun(@times,x_gcp,f_c.*f_Nc); X_gcp = bsxfun(@times,X_gcp,f_Nc);
  x_gcp = bsxfun(@times,x_gcp,f_p.*f_Np); X_gcp = bsxfun(@times,X_gcp,f_Np);
  
  % PROJECTION
  
  fprintf('Calculating p-value using 2D Projection method...  ');
  
  null_score_boost = 3;
  min_effect_size = 1.25;
  convolution_numbins = 1000;
  
  G.p = nan(ng,1);
  
  % Fan
  P1_gp = zeros(np, ng);
  
  for g=1:ng, if ~mod(g,1000), fprintf('%d/%d ',g,ng); end
    
    % STEP 1
    % for each sample, prioritize mutation categories according to how likely
    % it would be for this gene x sample to have a mutation there by chance.
    
    N = reshape(N_nonsilent(g,1:ncat,:),ncat,np)';
    n = reshape(n_nonsilent(g,1:ncat,:),ncat,np)';
    x = reshape(x_gcp(g,1:ncat,:),ncat,np)';
    X = reshape(X_gcp(g,1:ncat,:),ncat,np)';
    P0 = hyge2pdf(0,N,x,X);
    P1 = hyge2pdf(1,N,x,X);
    
    P1_gp(:, g) = 1 - P0(:,1) .* P0(:,2);
    
    % determine each patient's priority order of categories (according to P1)
    % left column of "priority" = least extreme category of mutation
    % right column of "priority" = most extreme category of mutation
    [tmp priority] = sort(P1,2,'descend');
    % sort the P arrays to put the columns in least->most priority order
    shft = (priority - repmat(1:ncat,np,1));
    map = reshape(1:(np*ncat),np,ncat);
    newmap = map + shft*np;
    P0 = P0(newmap);
    P1 = P1(newmap);
    P2 = 1-(P0+P1);  % note, P2 means "P(2+)"
    P2(P2<0) = 0;
    
    % STEP 2
    % for each sample, compute probability that it would have been of each (2-dimensional) degree.
    % degree=(d1,d2), where d=0 (no mut) ..... ncat (most extreme mut)
    % d1 is the MOST extreme mutation (or no mutation)
    % d2 is the SECOND MOST extreme mutation (or no mutation)
    % d1 can be 0-ncat; d2 can be 0-d1
    
    Pdeg = zeros(np,ncat+1,ncat+1);
    for d1=0:ncat, for d2=0:d1
      % has to have 0 in any/all categories > d1
      p = prod(P0(:,d1+1:end),2);
      if (d1>0)  % and (if d1>0)
        if (d1==d2)
          % if d1==d2, has to have 2+ in category d1
          p = p .* P2(:,d1);
        else
          % else:   has to have exactly 1 in category d1
          %         has to be clear in any/all categories (d2+1) to (d1-1)
          %         and (if d2>0) have (1 or 2+) in category d2
          p = p .* P1(:,d1);
          p = p .* prod(P0(:,d2+1:d1-1),2);
          if (d2>0)
            p = p .* (P1(:,d2)+P2(:,d2));
          end
        end
      end
      Pdeg(:,d1+1,d2+1) = p;
    end,end

    %% STEP 2a: calculate score for a sample being of each possible degree
    %% (uses new style, where score = -log10 probability of having a mutation in that category
    %% (zero score for no mutations)
    Sdeg = zeros(np,ncat+1,ncat+1);
    for d1=1:ncat, for d2=0:d1
      if d1==d2
        p = P2(:,d1);
      else
        if d2>0
          p = P1(:,d1).*P1(:,d2);
        else
          p = P1(:,d1);
        end
      end
      Sdeg(:,d1+1,d2+1) = -log10(p);
    end,end

    % null score boost
    priority2 = [zeros(np,1) priority];
    Sdeg(priority2==null_categ) = Sdeg(priority2==null_categ) + null_score_boost;
    
    % STEP 3
    % determine actual (two-dimensional) degree and score for each sample
    % sum scores to get score_obs for gene
    
    degree = zeros(np,2);
    score_obs = 0;
    for p = 1:np
      i = 1;
      for d = ncat:-1:1
        c = priority(p,d);
        if i==1
          if n(p,c)>=2
            degree(p,:) = [d d];
            i = 3;
          elseif n(p,c)==1
            degree(p,i) = d;
            i=i+1;
          end
        elseif i==2
          if n(p,c)>=1
            degree(p,i) = d;
            i=i+1;
          end
        else % i>2: done
          break
        end
      end
      score_sample = Sdeg(p,degree(p,1)+1,degree(p,2)+1);
      score_obs = score_obs + score_sample;
    end

    % minimum effect size
    score_obs = score_obs / min_effect_size;
    
    % for zero score, don't bother doing convolutions
    if score_obs<=0, G.p(g)=1; continue; end
    
    % STEP 4
    % compute P value for gene by convolutions
    
    numbins = convolution_numbins;
    binsize = score_obs / numbins;
    H = zeros(numbins,1);
    H(1) = 1;  % initial condition: all probability is in first bin
    
    % sequential convolution
    offset = min(numbins, round(Sdeg/binsize));
    ncols = (ncat+1)*(ncat+2)/2;
    newH = zeros(numbins,ncols);
    for p=1:np
      newH(:) = 0;
      col=1;
      for d1=0:ncat, for d2=0:d1
        o = offset(p,d1+1,d2+1);
        newH(o+1:end,col) = Pdeg(p,d1+1,d2+1) .* H(1:end-o);
        col=col+1;
      end,end
      H = sum(newH,2);
    end
    
    % save p-value
    G.p(g) = max(0,1-sum(H));
    
  end, fprintf('\n');   % next gene
  
  % Fan
  fid = fopen('all_genes.txt', 'wt');
  if fid>0
      for k=1:size(G.gene,1)
          fprintf(fid, '%s\n', G.gene{k, :});
      end
      fclose(fid);
  end
  dlmwrite(strrep(output_file, '.sig_genes.txt', '.p1gp'), P1_gp, 'delimiter', '\t', 'precision', '%.6f');
  
  % FDR
  G.q = calc_fdr_value(G.p);
  
  
  G = sort_struct(G,'p');
  save_struct(G,output_file);
  
  fprintf('Done.  Wrote results to %s\n',output_file);
  
end
% end of MutSig_runCV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% math subfunctions                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p=hyge2pdf(k,n,k1,n1)
  p = exp(gammaln(n1+2) - gammaln(k1+1) - gammaln(n1-k1+1) + ...
          gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ...
          gammaln(k1+k+1) + gammaln(n+n1-k-k1+1) - gammaln(n+n1+2));
  p = max(0,min(1,p));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p=hyge2cdf(k,n,k1,n1)
  p=0;
  for ki=0:k, p=p+hyge2pdf(ki,n,k1,n1); end
  p = max(0,min(1,p));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fdr=calc_fdr_value(p)
  if isempty(p)
    fdr = p;
    return
  end
  if size(p,1)==1
    p=p';
    trans=1;
  else
    trans=0;
  end
  [sp,ord]=sort(p);
  fdr=sp*size(p,1)./repmat((1:size(p,1))',1,size(p,2));
  fdr(fdr>1)=1;
  fdr=[fdr; ones(1,size(fdr,2))];
  for i=size(p,1):-1:1
    fdr(i,:)=min(fdr(i:(i+1),:),[],1);
  end
  fdr=fdr(1:(end-1),:);
  ordmat=ord+repmat(0:size(p,1):size(p,1)*(size(fdr,2)-1),size(ord,1),1);
  fdr(ordmat(:))=fdr(:);
  if trans
    fdr=fdr';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = hist3d(a,b,c,firsta,lasta,firstb,lastb,firstc,lastc)
  if nargin~=9, error('requires 9 input arguments'); end
  h = zeros(lasta-firsta+1,lastb-firstb+1,lastc-firstc+1);
  for bi=firstb:lastb, for ci=firstc:lastc
      h(:,bi,ci) = histc(a(b==bi & c==ci),firsta:lasta);
  end,end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% technical subfunctions: file input/output, etc.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ok = demand_file(fname)
  if ~iscell(fname), fname = {fname}; end
  ok = true(size(fname));
  for i=1:numel(fname), if ~mod(i,100), fprintf('%d/%d ',i,numel(fname)); end
    if ~exist(fname{i},'file'), ok(i)=false; end
  end, if numel(fname)>=100, fprintf('\n'); end
  nnf = sum(~ok);
  nf = sum(ok);
  if nargout==0
    if nf>=10, fprintf('%d files found\n',nf); end
  end
  if nnf>0
    blank = 0;
    for i=1:numel(fname)
      if ~ok(i)
        if strcmp(fname{i},'')
          blank=blank+1;
        else
          if nnf<20
            fprintf('\tNot found: %s\n',fname{i});
          end
        end
      end
    end
    if blank>0, fprintf('\tNot found: <blank> (%d)\n',blank); end
    if nargout==0
      if nnf==1, error('1 file not found'); else error('%d files not found',nnf); end
    end
  end
  if nargout==0, clear ok; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=dlmsep(s,d)
  if nargin==1
    d=9; % tab
  end

  pos=find(ismember(s,d));
  if ~isempty(pos)
    pos=[ 0 pos length(s)+1];
    for i=1:(length(pos)-1)
      res{i}=s((pos(i)+1):(pos(i+1)-1));
    end
  else
    res{1}=s;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = ensure_fopen(varargin)
  result = fopen(varargin{:});
  if (result==-1), error('ERROR WRITING TO %s',varargin{1}); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d=find_dlm(fname,dlm)
  if ~exist('dlm','var') || isempty(dlm)
    dlm=[ char(9) ',|'];
  end
  f=fopen(fname,'r');
  l=fgetl(f);
  for i=1:length(dlm)
    h(i)=length(find(l==dlm(i)));
  end
  [hm,hi]=max(h);
  d=dlm(hi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = flip(A)
  if sum(size(A)>1)>1, error('flip only works for vectors'); end
  if size(A,1)==1, A=fliplr(A); else A=flipud(A); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function numcols = get_colcount(fname,num_header_lines)
  if ~exist('num_header_lines','var'), num_header_lines = 1; end
  if ~exist(fname,'file'), error('%s not found',fname); end
  
  f = fopen(fname);
  for i=1:num_header_lines+1; l = fgetl(f); end
  numcols = sum(l==char(9))+1;
  fclose(f);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x=as_column(x)
  if size(x,2)>1
    x=x';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = impose_default_value(S,field,value,acceptable_values)
  if ~isfield(S,field) || isempty(getfield(S,field))
    if ischar(value) && strcmp(value,'*required*')
      error('%s is a required field of P',field);
    else
      try
        S=setfield(S,field,value);
      catch me
        fprintf('Error setting field "%s"\n',field);
        disp(me);disp(me.message);
      end
    end
  end
  
  if exist('acceptable_values','var')
    av = acceptable_values;
    v = getfield(S,field);
    if ischar(v) && isnumeric(av)
      v = str2double(v);
      S = setfield(S,field,v);
    end
    if ~ischar(v) && ischar(av)
      error('%s is assigned value of different type than acceptable_values',field);
    end
    try
      ism = ismember(v,av);
    catch me
      error('%s is assigned value of different type than acceptable_values',field);
    end
    if ~ism
      fprintf('Acceptable values for %s:\n',field); disp(av);
      fprintf('Attempted to set to:\n'); disp(v)
      error('Invalid setting of %s',field);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m = listmap(a,b)
  if ischar(a), a={a}; end
  if ischar(b), b={b}; end
  m = nan(length(a),1);
  [a ai aj] = unique(a);
  [c ia ib] = intersect(a,b);
  for i=1:length(c), if ~mod(i,1e5), fprintf('%d/%d ',i,length(c)); end
    m(aj==ia(i)) = ib(i);
  end, if length(c)>=1e5, fprintf('\n'); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = load_lines(fname)
  X = load_textfile(fname);
  L = text_to_lines(X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = load_struct(varargin)
% load_struct(filename, format, separator_character, header_lines, <other_parameters>, P)
  
% loads a tab-delimited table from the specified file
% into a structure with field names based on the column headings
%
% piggybacks on read_table()
% all parameters are passed to read_table(),
% except last parameter if it is a struct (i.e. P struct)
% P struct can have following members:
%   P.lowercase_fieldnames: if true or 1, fieldnames are converted to lowercase
%
% WAS: load_struct(filename, format, header_lines, lowercase_fieldnames)
%   format is [] by default.
%   header_lines is 1 by default.
%   if header_lines is set to 0, field names are col1, col2, etc.
%   if lowercase_fieldnames is 1, fieldnames are converted to lowercase
%
% Mike Lawrence 2008-04-29
% fixed 2010-12-10 to give all its parameters to read_table
%   (except last parameter if a P struct)

% check parameters

  args = varargin;
  if length(args)<1
    error('filename required');
  end
  filename = args{1};
  if ~ischar(filename)
    error('first parameter should be filename (character string)');
  end
  demand_file(filename);
  
  % see if last parameter is a P struct
  if isstruct(args{end})
    P = args{end};
    args = args(1:end-1);
  else
    P = [];
  end
  
  P = impose_default_value(P,'lowercase_fieldnames',false);
  P = impose_default_value(P,'ignore_poundsign_lines',true);
  
  if length(args)>=2 
    % ok
  end
  if length(args)>=3
    if isempty(args{3})
      error('third parameter should be separator character, e.g. char(9)');
    end
    if isnumeric(args{3})
      error('header_lines is now the fourth parameter... please update call to load_struct');
    end
  end
  if length(args)>=4
    if islogical(args{4})
      error('lowercase_fieldnames has been moved to P struct... please update call to load_struct');
    end
    if ~isnumeric(args{4})
      error('fourth parameter should be number of header lines');
    end
  end
  
  % HANDLE COMMENT LINES:
  default_header_lines = 1;
  %% see if table has comment lines at the beginning (start with #): if so, increment header_lines to skip them
  n_comment_lines = 0;
  
  if P.ignore_poundsign_lines
    f = fopen(args{1});
    while(true)
      x = fgetl(f);	
      if (isempty(x))||(x(1)=='#') % skip empty lines or lines starting with #
        n_comment_lines = n_comment_lines + 1;
        continue
      elseif strncmp(x,'Oncotator v',11)
        fprintf('Un-poundsigned Oncotator header detected and skipped.\n');
        n_comment_lines = n_comment_lines + 1;
        continue
      else
        break
      end
    end
    fclose(f);
  end
  
  default_header_lines = 1 + n_comment_lines;
  
  %% default args
  if length(args)==1, args = [args {''}]; end         % format string
  if length(args)==2, args = [args {char(9)}]; end       % separator character
  if length(args)==3, args = [args {default_header_lines}]; end          % number of header lines
  
  % default whitespace
  has_already = false;
  for i=4:length(args)
    if ischar(args{i}) && strcmpi(args{i},'whitespace'), has_already=true; break; end
  end
  if ~has_already, args = [args {'whitespace'} {'\b\r'}]; end
  
  % default bufSize
  if verLessThan('matlab','8.4')
    has_already = false;
    for i=4:length(args)
      if ischar(args{i}) && strcmpi(args{i},'bufSize'), has_already=true; break; end
    end
    if ~has_already, args = [args {'bufSize'} {50000}]; end
  end
  
  % LOAD TABLE
  try
    table = read_table(args{:});
    nf = length(table.dat);
  catch me
    q = load_lines(args{1});
    if isempty(q)
      fprintf('\n%s is a blank file\n',args{1});
      table = [];
      table.dlm  = args{3};
      table.headers = {{}};
      table.dat = {{}};
      nf = 0;
    else
      disp(me);
      disp(me.message);
      error('Error loading struct file');
    end
  end
  
  if isempty(table.headers)
    table.headers{1} = cell(nf,1);
    for f=1:nf
      table.headers{1}{f} = sprintf('col%d', f);
    end
  end
  
  % process header line
  
  fields = table.headers{end};
  if length(fields)~=nf
    fprintf('Header line has %d column headers instead of the expected %d:\n',length(fields),nf);
    fields{:}
    error('Unable to parse table header line.');
  end
  
  % remove illegal characters from column headings
  % and convert to list of unique field names
  
  if P.lowercase_fieldnames, fields = lower(fields); end
  fields = regexprep(fields, '\W','');   % remove any characters except A-Z, a-z, 0-9, underscore
  fields_orig = fields;
  fields = genvarname(fields_orig);
  
  % preserve "end", because it's only going to be a field name, not a variable name
  for f=1:nf
    if strcmp(fields_orig{f}, 'end')
      fields{f} = 'end';
      break
    end
    if strcmp(fields_orig{f}, 'End')
      if P.lowercase_fieldnames
        fields{f} = 'end';
      else
        fields{f} = 'End';
      end
      break
    end
  end
  
  S = struct();
  for f=1:nf
    S = setfield(S, fields{f}, table.dat{f});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = load_struct_specify_string_cols(fname,string_cols,num_header_lines,lowercase_fieldnames)

  if ~exist('fname','var'), error('Must specify fname'); end
  if ~exist('string_cols','var'), string_cols = []; end
  if ~exist('num_header_lines','var'), num_header_lines = 1; end
  if ~exist('lowercase_fieldnames','var'), lowercase_fieldnames = false; end
  
  if ~exist(fname,'file'), error('%s not found',fname); end
  
  numcols = get_colcount(fname,num_header_lines);
  
  is_string = false(1,numcols);
  is_string(string_cols) = true;
  format = [];
  for i=1:numcols
    if is_string(i), format = [format '%s'];
    else format = [format '%f']; end
  end
  
  P=[]; P.lowercase_fieldnames = lowercase_fieldnames;
  X = load_struct(fname,format,char(9),num_header_lines,P);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t = load_textfile(filename)
  in = fopen(filename);
  t = fread(in,'uint8=>char')';
  fclose(in);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = make_boolean(S,varargin)
  fields = {};
  for i=1:length(varargin)
    fields = [fields varargin{i}];
  end
  for i=1:length(fields)
    if isfield(S,fields{i})    
      x = getfield(S,fields{i});
      if ~islogical(x) && ~isnumeric(x)
        x = str2double(x);
      end
      if ~islogical(x)
        x = (x~=0);
      end
      S = setfield(S,fields{i},x);
    else
      fprintf('No such field: %s\n', fields{i});
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = make_numeric(S,varargin)
  fields = {};
  for i=1:length(varargin)
    fields = [fields varargin{i}];
  end
  for i=1:length(fields)
    if isfield(S,fields{i})    
      x = getfield(S,fields{i});
      if isnumeric(x) || islogical(x)
      else
        x = str2double(x);
        S = setfield(S,fields{i},x);
      end
    else
      fprintf('No such field: %s\n', fields{i});
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = nansub(X,idx,filler)

  if numel(X)==2 && size(X,1)==1 && size(X,2)>1
    %   fprintf('note: converting first argument to column vector\n');
    X = X';
  end

  if iscellstr(X) && size(X,1)==1 && size(X,2)>1
    X=X';
  end

  if islogical(X)
    type = 0;
  elseif isnumeric(X)
    type = 1;
  elseif iscell(X)
    type = 2;
  else
    error('Unsuuported array type');
  end
  
  if ~exist('filler','var')
    if type==0
      filler = false;
    elseif type==1
      filler = nan;
    elseif type==2
      filler = {''};
    else
      error('Inconsistent behavior with "type"');
    end
  end
  
  if type==0
    if ~islogical(filler)
      error('Inappropriate filler for logical array');
    end
  elseif type==1
    if ~isnumeric(filler)
      error('Inappropriate filler for numeric array');
    end
  elseif type==2
    if ischar(filler)
      filler = {filler};
    end
    if ~iscell(filler)
      error('Inappropriate filler for cell array');
    end
  else
    error('Inconsistent behavior with "type"');
  end
  
  sz = size(X); sz(1) = length(idx);
  Y = repmat(filler,sz);
  idx2 = find(~isnan(idx) & idx>=1 & idx<=length(X));
  Y(idx2,:,:,:,:) = X(idx(idx2),:,:,:,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function st=newline
  if ispc
    st='\r\n';
  else
    st='\n';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = parse(S,r,f,numeric)
  if ~iscell(S), S = tolines(S); end
  tokens = regexp(S,r,'tokens');
  if ~iscell(f), f = {f}; end
  x = tokens2struct(tokens,f);
  if exist('numeric','var'), x = make_numeric(x,f(numeric)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fl,fid]=read_dlm_file(fname,dlm,nlines)
% READ_DLM_FILE  Read all or some of a delimited file into cell arrays-of-arrays.
% [FL,FID] = READ_DLM_FILE(FNAME,DLM,NLINES)
%    Reads a file (FNAME) and partitions to line using the delimeter (DLM).
%    The default DLM value is \t (9). The optional NLINES parameter
%    limits the number of lines read. 
%    
%    FL is a cell array of cell arrays of strings.  FL{n} is a cell array 
%    strings containing the nth line of the file.  
%    FID is the fid for the opened file, FNAME.
%
%
% Gaddy Getz
% Cancer Genomics
% The Broad Institute
% gadgetz@broad.mit.edu
%

  if nargin==1
    dlm=9;
  end
  
  if ischar(fname)
    fid=fopen(fname,'r');
  else
    fid=fname;
  end
  
  ln=1;
  if ~exist('nlines','var')
    nlines=Inf;
    do_close=1;
  else
    do_close=0;
  end
  
  had_output=0;
  fl={};
  while(ln<=nlines)
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    fl{ln}=dlmsep(tline,dlm);         
    ln=ln+1;
    if mod(ln,1000)==0
      verbose(['...' num2str(ln)],30);
      had_output=1;
    end
  end
  if do_close
    fclose(fid);
    fid=-1;
  end
  ln=ln-1;
  if had_output
    verbose([newline],30);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tab=read_table(fname,format,dlm,headerlines,varargin)

  fpos=0;
  if headerlines~=0
    if ischar(fname)
      f=fopen(fname,'r');
    else
      f=fname;
      fpos=ftell(f);
    end
    if length(dlm)~=1
      tab.dlm=find_dlm(fname,dlm);
    else
      tab.dlm=dlm;
    end
    if headerlines>0
      tab.headers=read_dlm_file(f,dlm,headerlines);
    elseif headerlines==-1 % R convention
      headerlines=1;
      tab.headers=read_dlm_file(f,dlm,headerlines);
      tab.headers{1}=['EMPTY' tab.headers{1,:}];
    end
  else
    if ischar(fname)
      f=fopen(fname,'r');
    else
      f=fname;
      fpos=ftell(f);
    end
    tab.headers={};
    tab.dlm=dlm;
  end
  
  if isempty(format)
    if isempty(tab.headers)
      error('must have either format or headerlines');
    else
      format=[repmat('%s',1,length(tab.headers{end})) '\n'];   % (to allow for multiple header lines)
    end
  elseif iscell(format)
    if isempty(tab.headers)
      error('must have either format or headerlines');
    else
      if length(format)==1
        format=[repmat(format{1},1,length(tab.headers{end})) '\n'];
      else
        format=[format{1} repmat(format{3},1,length(tab.headers{end})-format{2}) '\n'];
      end
    end
  end

  if strcmp(format((end-1):end),'\n')
    format=format(1:(end-2));
  end
  
  verbose(['Reading file using format:' format],10);
  fseek(f,fpos,'bof');
  tab.dat=textscan(f,format,'headerLines',headerlines,'delimiter',tab.dlm,varargin{:});
  fclose(f);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s order]=reorder_struct(s,order)
  
  if nargin~=2, error('reorder_struct(s,order)'); end
  
  if islogical(order), order = find(order); end
  if ischar(order)
    if strcmpi(order,'end')
      order = slength(s);
    else
      error('invalid index parameter');
    end
  end
  
  order = as_column(order);
  nanflag = any(isnan(order));
  fields = fieldnames(s);
  nf = length(fields);

  for i=1:nf
    f = getfield(s,fields{i});
    if nanflag
      f = nansub(f,order);
    else
      f = f(order,:,:,:,:,:,:,:,:,:);
    end
    s = setfield(s,fields{i},f);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = save_struct(S,filename,noheader,deprecated_formatflag)
%
% save_struct(S, filename)
%
% writes a tab-delimited table from the given struct.
% writes a header line based on the field names.
%
% fields are written in the order in which they appear in the struct.

  if nargin==0, error('requires an input argument'); end
  if nargin>4, error('too many input arguments'); end
  if nargin==4, fprintf('USE OF DEPRECATED FORMATFLAG\n'); end
  
  if nargout==0 && nargin<2, error('requires at least 2 arguments'); end
  if nargout==1 && (nargin<2 || isempty(filename))
    return_string_only = true;
  else
    return_string_only = false;
  end
  if nargout>1, error('only a single output possible'); end
  
  if exist('noheader','var') && ~isempty(noheader) && (~islogical(noheader) || ~(noheader==false))
    if (islogical(noheader) && noheader==true) || strncmpi(noheader,'no_header',9) || strncmpi(noheader,'noheader',8)
      noheader = true;
    else
      error('third parameter should be "no_headers", true, false, empty, or nothing.');
    end
  else
    noheader = false;
  end
  
  if ~return_string_only
    out = ensure_fopen(filename,'wt');
  end
  
  if ~isstruct(S)
    if isempty(S)
      S = struct;
    else
      error('S should be a struct');
    end
  end
  
  fld = fieldnames(S);
  nf = length(fld);
  slen = slength(S);
  
  % see if struct is empty
  if slen==0
    if return_string_only
      F = '';
    else
      for i=1:nf
        fprintf(out,fld{i});
        if i==nf
          fprintf(out,'\n');
        else
          fprintf(out,'\t');
        end
      end
      fclose(out);
    end
    return
  end

  % see if struct is to long too handle all at once
  chunksize = round(1e7/nf);
  if slen>chunksize
    if return_string_only
      error('struct is too large to convert all at once in memory');
    end
    
    for st=1:chunksize:slen
      fprintf('HUGE STRUCT: SAVING CHUNK %d/%d\n', ceil(st/chunksize), ceil(slen/chunksize));
      en = min(slen,st+chunksize-1);
      Si = reorder_struct(S,st:en);
      if exist('deprecated_formatflag','var')
        F = save_struct(Si,[],noheader,deprecated_formatflag);
      else
        F = save_struct(Si,[],noheader);
      end
      fwrite(out,F);
      clear F;
      noheader = 'no_header'; % subsequent chunks omit header
    end
    
  else   % struct is not too big to save all at once
    F = cell(1,nf*2);
    nr = -1;
    tt0 = tic;
      
      for f=1:nf
        tt1 = toc(tt0);
        if tt1>10
          if ~exist('flag01','var')
            fprintf('  [save_struct] ');
            flag01 = true;
          end
          fprintf('%d/%d ',f,nf);
        end
        C = getfield(S,fld{f});    % get column
                                   % check for legal type
        if isempty(C)
          C = {};
        elseif isnumeric(C) || islogical(C)
          if ndims(C)>2
            fprintf('Field %s is multidimensional: skipping\n', fld{f});
            C = {};
          elseif size(C,2)>1
            fprintf('Field %s is not a column vector: skipping\n', fld{f});
            C = {};
          else
            if exist('deprecated_formatflag','var') && deprecated_formatflag     % convert to cell array of strings
              C = cellstr(num2str(C));
            else                       % note: "-" is important to avoid extra spaces in output
              if any(mod(C,1))
                C = cellstr(num2str(C,'%-d'));     % not all integers
              else
                C = cellstr(num2str(C,'%-.0f'));    % all integers
              end
            end
          end
        else  % column is a cell
          idx = find(cellfun(@iscell,C));
          if ~isempty(idx)
            fprintf('WARNING: Field %s contains %d entries that are cells:\n', fld{f}, length(idx));
            fprintf('         replacing these with "?????{cell}?????"\n');
            C(idx) = repmat({'?????{cell}?????'},length(idx),1);
          end
          idx = find(~cellfun(@isempty,C) & ~cellfun(@ischar,C));
          if ~isempty(idx)
            fprintf('WARNING: Field %s contains %d entries that are not chars:\n', fld{f}, length(idx));
            fprintf('         replacing these with "?????{unknown_format}?????"\n');
            C(idx) = repmat({'?????{unknown_format}?????'},length(idx),1);
          end
        end
        if isempty(C)
          if nr==-1, error('Problematic leftmost column'); end
          C = repmat({'...'},nr,1);
        end
        % check for consistent column length
        if nr==-1, nr=length(C); end
        if nr~=length(C), error('Field %s is a different length', fld{f}); end
        % add column title
        if ~noheader, C = [fld{f}; C]; end
        % add column to file
        F{f*2-1} = C;
        % add tab or newline
        F{f*2} = repmat({char(9+(f==nf))},length(F{f*2-1}),1);
      end
      
      if tt1>10, fprintf(' [collapse]'); end
      F = strcat(F{:});              % collapse columns to lines
      F = [F{:}];                    % collapse lines to file
      
      if ~return_string_only
        if tt1>10, fprintf(' [write]'); end
        fwrite(out,F);
      end
      
      if tt1>10, fprintf('\n'); end
  end
  
  if ~return_string_only
    fclose(out);
    clear F
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = rename_field(S, oldname, newname)

  if iscell(oldname) && iscell(newname)
    if ~iscell(newname) || length(oldname)~=length(newname), error('lists must be same length'); end
  elseif ~iscell(oldname) && ~iscell(newname)
    oldname = {oldname};
    newname = {newname};
  else
    error('improper parameters');
  end
  
  flds = fieldnames(S);
  
  for i=1:length(oldname)
    f = getfield(S, oldname{i});
    S = setfield(S, newname{i}, f);
    if ~strcmp(oldname{i},newname{i})
      S = rmfield(S, oldname{i});
    end
    idx = find(strcmp(flds,oldname{i}));
    if length(idx)~=1, error('unexpected behavior'); end
    flds{idx} = newname{i};
  end
  
  S = order_fields_first(S,unique_keepord(flds));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = order_fields_first(S,first_flds)

  if ischar(first_flds), first_flds = {first_flds}; end
  
  all_flds = fieldnames(S);
  
  if ~isempty(setdiff(first_flds,all_flds)), error('Some of those fields don''t exist'); end
  
  rest_flds = all_flds;
  rest_flds(ismember(rest_flds,first_flds)) = [];
  
  S = orderfields(S,[as_column(first_flds);as_column(rest_flds)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u ui uj] = unique_keepord(x,varargin);
  if exist('varargin','var') && length(varargin)>=1 && ischar(varargin{1}) && (strcmpi(varargin{1},'first')|strcmpi(varargin{1},'last'))
    error('please do not specify "first" or "last" with this function.  (default is "first")');
  end
  
  [u1 ui1 uj1] = unique(x,'first',varargin{:});
  
  [ui ord] = sort(ui1);
  u = x(ui1(ord));
  [tmp ord2] = sort(ord);
  uj = ord2(uj1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function l = slength(S)
  l=NaN;
  if isstruct(S)
    l = 0;
    if ~isempty(S) && ~isempty(fieldnames(S))
      f = fields(S);
      nf = length(f);
      len = nan(nf,1);
      for i=1:nf
        f1 = getfield(S,f{i});
        if ischar(f1), len(i) = 1;
        else len(i) = size(f1,1); end
      end
      ulen = unique(len);
      if length(ulen)==1, l = ulen;
      else
        fprintf('Warning: deprecated use of slength for structure with fields of nonuniform length\n');
        l = len(1);
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s2,ord]=sort_struct(s1,keyfield,order)
  if length(keyfield)==0, return; end
  if ~iscell(keyfield)
    keyfield = {keyfield};
  end
  if ~exist('order','var')
    order = repmat(1,length(keyfield),1);
  end
  if length(order) ~= length(keyfield)
    error('order and keyfield must have same number of elements');
  end
  
  orig_len = slength(s1);
  s2=s1;
  ord=(1:orig_len)';
  fields = fieldnames(s1);
  nf = length(fields);
  
  for k=length(keyfield):-1:1
    f = getfield(s2,keyfield{k});
    if length(f)<orig_len, error('Attempted to sort on truncated field "%s"',keyfield{k}); end
    if order(k)==1
      if isnumeric(f)
        [tmp ordi] = sortrows(f);
      else
        [tmp ordi] = sort(f);
      end
    elseif order(k)==-1
      if isnumeric(f)
        [tmp ordi] = sortrows(f,-1);
      else
        [tmp ordi] = sort(flip(f));
      end
    else
      error('Unknown order %d',order(k));
    end
    for i=1:nf
      f = getfield(s2,fields{i});
      f = f(ordi,:,:,:,:,:,:,:,:,:);
      s2 = setfield(s2,fields{i},f);
    end
    ord = ord(ordi);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = split(string, delim)
  if ischar(string), string = {string}; charflag=true; else charflag=false; end
  ns = length(string);
  R = cell(ns,1);
  for z=1:ns, if mod(z,10000)==0, fprintf('%d/%d ',z,ns); end
    
    dpos = [0 find(string{z}==delim) length(string{z})+1];
    nt = length(dpos)-1;
    R{z} = cell(nt,1);
    for t=1:nt
      R{z}{t} = string{z}(dpos(t)+1:dpos(t+1)-1);
    end
  end, if ns>=10000, fprintf('\n'); end
  if charflag, R = R{1}; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = text_to_lines(X)
  while ~isempty(X)
    if X(end)==10, X(end)=[];   % remove trailing blank lines
    else break; end
  end
  if isempty(X), L={};
  else L = split(X,char(10)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = tokens2struct(T,fieldnames)
  nt = length(T);
  nf = length(fieldnames);
  S = [];
  for f=1:nf
    F = cell(nt,1);
    for t=1:nt
      if isempty(T{t}), F{t} = '';
      else F{t} = T{t}{1}{f}; end
    end
    S = setfield(S,fieldnames{f},F);
  end
end

function a = tolines(a)
  a = split(a,char(10));
  a = a(~cellfun('isempty',a));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function verbose(str,level,varargin)

  global VERBOSE_LEVEL
  global VERBOSE_FILE

  if nargin==1
    level=1;
  end
  
  if isempty(varargin)
    
    str = char(regexprep(cellstr(str),'%','%%'));
    if ~isempty(VERBOSE_LEVEL) && (level<=VERBOSE_LEVEL)
      fprintf(1,[str repmat('\n',size(str,1),1)]');  %escape the % to prevent line from commenting
      if ~isempty(VERBOSE_FILE)
        if ~exist(VERBOSE_FILE,'file')
          fid = fopen(VERBOSE_FILE,'w');
        else
          fid = fopen(VERBOSE_FILE,'a');
        end
        
        fprintf(fid,str,varargin{:});  %escape the % to prevent line from commenting
        fprintf(fid,'\n');
        fclose(fid);
      end
    end
    
  else
    
    if ~isempty(VERBOSE_LEVEL) && (level<=VERBOSE_LEVEL)
      fprintf(1,str,varargin{:});  %escape the % to prevent line from commenting
      fprintf(1,'\n')
      if ~isempty(VERBOSE_FILE)
        fid = fopen(VERBOSE_FILE,'a');
        fprintf(fid,str,varargin{:});  %escape the % to prevent line from commenting
        fprintf(fid,'\n');
        fclose(fid);
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function c = assign_65x4_to_categ_set(K)
% c = assign_65x4_to_categ_set(K)
%
% Given a set of k categories as struct K with the following fields:
%   left   = subset of 'ACGT', representing 5' base
%   right  = subset of 'ACGT', representing 3' base
%   from   = subset of 'AC', representing mutated base (after strand collapse)
%   change = subset of 'tfs', representing Transition, Flip transversion, Skew transversion
%   type   = either "point" or "non-point"
%
% LEGACY CASE:  if name or type contains "null" or "indel" (or "double_null"), then type is taken to be "non-point"
%
% Maps the set of 65 territory categories and 4 newbases
%    onto this reduced set of categories.
%
% Returns matrix c, with 65 rows, k columns, and 4 pages (one for each newbase);
%   cells with "1" are counted toward that category.
%   For indel/null categories, all pages (newbases) are set to the same value.
%
% Mike Lawrence 2010-01-27

require_fields(K,{'left','right','from','change','type'});
nk = slength(K);
% legacy cases
idx = grepi('indel|null',K.type,1);
if ~isempty(idx), K.type(idx) = repmat({'non-point'},length(idx),1); end
if isfield(K,'name')
  idx = grepi('indel|null',K.name,1);
  if ~isempty(idx), K.type(idx) = repmat({'non-point'},length(idx),1); end
end

X = generate_categ_context65_names();
require_fields(X,{'num','name'});
X = sort_struct(X,'num');  % (already sorted, but doesn't hurt to make sure)
X = parse(X.name,'(.) in (.)_(.)',{'from','left','right'});

base = 'ACGT';
complement(base) = 'TGCA';
whatchange('A','ACGT') = 'nstf';
whatchange('C','ACGT') = 'snft';
whatchange('G','ACGT') = 'tfns';
whatchange('T','ACGT') = 'ftsn';

c = zeros(65,nk,4);

for x=1:65   % for each context65 
  from = X.from{x};
  if isempty(from), continue; end   % "N"
  left = X.left{x};
  right = X.right{x};
  if ismember(from,'GT')
    from = complement(from);
    left = complement(X.right{x});
    right = complement(X.left{x});
  end
  for k=1:nk
    % decide if the context belongs to this category
    if ismember(from,K.from{k}) && ismember(left,K.left{k}) && ismember(right,K.right{k})
      % yes, it does
      if strcmpi(K.type{k},'non-point')
        c(x,k,:) = 1;   % for non-point mutations (indel/null), "newbase" doesn't matter
      elseif strcmpi(K.type{k},'point')
        for n=1:4     % which newbase
          oldbase = X.from{x};
          newbase = base(n);
          change = whatchange(oldbase,newbase);
          if ismember(change,K.change{k})
            c(x,k,n) = 1;
          end
        end
      else
        error('Unknown type: %s',K.type{k});
      end
    end     
  end  % next k
end  % next x

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function st=catn(s,dlm)

if ~exist('dlm','var')
  dlm=9;
end

if iscell(s)
  s=strvcat(s);
end
  
st=[num2str((1:size(s,1))') repmat(dlm,size(s,1),1) s];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = collapse_Nn_64_by_strand(in)
% rows:  64 categories in standard order
% columns:  [N ->A ->C ->G ->T]
% Mike Lawrence 2009-12-11

if size(in,1)~=64, error('input must have 64 rows'); end
if size(in,2)~=5, error('input must have 5 columns (N A C G T)'); end

out = in(1:32,:,:,:,:);
compbase('ACGT') = 'TGCA';
X = generate_categ_context65_names();
for i=1:32
  oldname = X.name{i};
  newname = [compbase(oldname(1)) ' in ' compbase(oldname(end)) '_' compbase(oldname(end-2))];
  j = find(strcmp(newname,X.name));
  out(i,1,:,:,:) = in(i,1,:,:,:) + in(j,1,:,:,:);
  out(i,2,:,:,:) = in(i,2,:,:,:) + in(j,5,:,:,:);
  out(i,3,:,:,:) = in(i,3,:,:,:) + in(j,4,:,:,:);
  out(i,4,:,:,:) = in(i,4,:,:,:) + in(j,3,:,:,:);
  out(i,5,:,:,:) = in(i,5,:,:,:) + in(j,2,:,:,:);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = collapse_Nn_65_to_32(in)
% rows:  65 categories in standard order
% columns:  [N ->A ->C ->G ->T]

if size(in,1)~=65, error('input must have 65 rows'); end
if size(in,2)~=5, error('input must have 5 columns (N A C G T)'); end

% (discard bottom row, "N")
out = collapse_Nn_64_by_strand(in(1:64,:));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = concat(strings,separator)
% concat(strings,separator)
%
% joins strings into one string, separated with the specified character/string

c='';
for i=1:length(strings)
  if ischar(strings(i))
    c=[c strings(i)];
  elseif isnumeric(strings(i))
    c=[c num2str(strings(i))];
  elseif isnumeric(strings{i})
    c=[c num2str(strings{i})];
  elseif iscell(strings(i))
    c=[c strings{i}];
  else
    error('concat does not know how to handle that kind of input');
  end
  if i<length(strings)
    c=[c separator];
  end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = concat_structs_keep_all_fields(X)

% get comprehensive list of fields
allflds = [];
type = [];  % 1 = numeric, 2 = boolean, 3 = cell
for i=1:length(X)
  fn = fieldnames(X{i});
  allflds = [allflds; fn];
  for j=1:length(fn)
    g = getfield(X{i},fn{j});
    if isnumeric(g), t = 1;
    elseif islogical(g), t=2;
    elseif iscell(g), t = 3;
    else error('Unknown type: operand %d field %s\n',i,fn{j});
    end
    type = [type; t];
  end
end
[flds ui uj] = unique(allflds,'first');

% make sure types are compatible
ot = type;
type = nan(length(flds),1);
for i=1:length(flds)
  t = ot(uj==i);
  if length(unique(type(i)))>1, error('Operands have "%s" of different types',flds{i}); end
  type(i) = t(1);
end

% for each operand, if it doesn't have the field in question, then add a blank version
for i=1:length(X)
  for j=1:length(flds)
    if ~isfield(X{i},flds{j})
      if type(j)==1, z = nan(slength(X{i}),1);
      elseif type(j)==2, z = false(slength(X{i}),1);
      elseif type(j)==3, z = repmat({''},slength(X{i}),1);
      else error('Inconsistent behavior');
      end
      X{i} = setfield(X{i},flds{j},z);
    end
  end
end

X = concat_structs(X);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function S = concat_structs(slist)
%
% concat_structs(slist)
%
% slist should be a cell array of structs,
%   all of which have the same fields
% 
% returns a new struct in which each field is a concatenation of
%   the corresponding fields from each of the input structs
%
% Mike Lawrence 2008-2010
%

if ~iscell(slist)
  error('input should be a cell array of structs');
end

if length(slist)==0
  S = {};
elseif length(slist)==1
  S = slist{1};
else

  ns = numel(slist);
  allflds = cell(ns,1);
  for i=1:ns
    if isempty(slist{i}), continue; end
    if ~isstruct(slist{i}), error('input should be a cell array of structs'); end
    allflds{i} = fieldnames(slist{i});
  end
  allflds = cat(1,allflds{:});
  [flds ai aj] = unique_keepord(allflds);
  h = histc(aj,1:length(flds));
  if ~all(h==ns)
    count(allflds);
    error('all input structs must have same fields.  use concat_structs_keep_all_fields.');
  end
  [tmp ord] = sort(ai);

  rflag = false;
  S = [];
  for fno=1:length(flds)
    type = nan(ns,1);
    f = cell(ns,1);
    for i=1:ns
      if isempty(slist{i}), continue; end
      f{i} = getfield(slist{i},flds{fno});
      if isnumeric(f{i}) || islogical(f{i}), type(i) = 1;
      elseif iscell(f{i}), type(i) = 2;
      else type(i) = 3;
      end
    end
    if any(type==3), error('Incompatible type encountered in %s',flds{fno}); end
    if any(type==1) && any(type==2)    
       if ~rflag, fprintf('Reconciling mixed cell+numeric fields for %s\n',flds{fno}); rflag=true; end
       for i=1:ns
         if type(i)==2
           try
             f{i} = str2double(f{i});
           catch me
             error('Unable to resolve mixed cell+numeric case encountered in %s',flds{fno});
    end,end,end,end
    S = setfield(S, flds{fno}, cat(1,f{:}));
  end

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [b,u] = count(a, sort_by_frequency, varargin)

if ~exist('sort_by_frequency', 'var'), sort_by_frequency = false; end

if nargin>=2 && length(sort_by_frequency)>1
  % probably meant to call xcount
  fprintf('Did you mean "xcount"?\n');
  xcount(a, sort_by_frequency, varargin{:})
  return
end

%% finds the unique elements of array and counts how many of each there

if size(a,1)==1
   a = a';
end

[u ui uj] = nanunique(a);
nu = length(u);

% make sure _u_ is cell array of strings

if ~iscell(a)
  tmp=u;
  u = cell(length(tmp),1);
  for i=1:length(tmp)
    if ischar(tmp(i))
      u(i) = {tmp(i)};
    else
      u(i) = {num2str(tmp(i))};
    end
  end
end

b = zeros(nu,1);
for j=1:nu
    b(j) = length(find(uj==j));
end

if sort_by_frequency == 1
  [b ord] = sort(b);
  u = u(ord);
elseif sort_by_frequency == -1
  [b ord] = sort(b, 'descend');
  u = u(ord);
end

bb = cell(nu,1);
for j=1:nu
    bb{j} = b(j);
end

u = [u; '----TOTAL'];
bb = [bb; {length(a)}];
nu=nu+1;

L=zeros(nu,1);
for i=1:nu
 L(i)=length(u{i});
end
maxL = max(L);

if nargout==0
  fprintf('\n');
  f = ['    %' num2str(maxL) 's: [%d]\n'];
  for i=1:nu
    fprintf(f, u{i}, bb{i});
  end
end

if nargout<2, clear u; else u=u(1:end-1); end
if nargout<1, clear b; end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t = direc(dirname)

if exist(dirname,'dir')
  if dirname(end)~='/', dirname = [dirname '/']; end
end

t = {};

tmp = rdir(dirname);
if isempty(tmp), return; end
t = list2cell(tmp.name);
t = remove_ansi(t);
t = grepv('^\.',t);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ensure_dir_exists(dirname)

if ~iscell(dirname)
  if ~exist(dirname,'dir'), mkdir(dirname); end
else
  for i=1:numel(dirname)
    if ~exist(dirname{i},'dir'), mkdir(dirname{i}); end
  end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ensure_writeable(fname)

try
  [path name ext] = fileparts(fname);
  if ~isempty(path), ensure_dir_exists(path); end
  testfile = [fname '.' rand() '.test'];
  save_textfile('test',testfile);
  delete(testfile);
catch me
  error('%s is not writeable',fname);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [b,a]=exchange_vars(a,b)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function names = find_good_names_for_mutation_categories(autonames)

if ~iscell(autonames)
  noncellflag = true;
  autonames = {autonames};
else
  noncellflag = false;
end

z.code = {'At','Af','As','Atf','Afs','Ats','Atfs',...
          'Ct','Cf','Cs','Ctf','Cfs','Cts','Ctfs',...
          'ACt','ACf','ACs','ACtf','ACfs','ACts','ACtfs'};
z.name = {'A->G','A->T','A->C','A->(G/T)','A->(T/C)','A->(C/G)','A->mut',...
          'C->T','C->G','C->A','C->(T/G)','C->(G/A)','C->(A/T)','C->mut',...
          'N->transit','N->flip','N->skew','N->nonskew','N->transver','N->nonflip','N->mut'};

p = parse(autonames,'^([ACGT]+)\[([AC])+->([tfs]+)\]([ACGT]+)$',{'before','at','change','after'});
p.muttype = mapacross(stringsplice([p.at p.change]),z.code,z.name);

for i=1:length(autonames)
  names{i,1} = regexprep(p.muttype{i},'^(.*)(->.*)$',[p.before{i} 'p*$1p' p.after{i} '$2']);
end

names = regexprep(names,'ACGTp','');
names = regexprep(names,'pACGT','');
names = regexprep(names,'^([ACGT])([ACGT])p','($1/$2)p');
names = regexprep(names,'^([ACGT])([ACGT])([ACGT])p','($1/$2/$3)p');
names = regexprep(names,'p([ACGT])([ACGT])->','p($1/$2)->');
names = regexprep(names,'p([ACGT])([ACGT])([ACGT])->','p($1/$2/$3)->');
names = regexprep(names,'^\*([ACN])->','$1->');
names = regexprep(names,'^N->','');

if ~noncellflag
  autonames = autonames{1};
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,Sc,Sv] = find_mut_categs(Nn,P)
% Nn should be a 32x5 table as from collapse_Nn_64_by_strand
%   rows:  32 strand-collapsed categories
%   columns:  [N ->A ->C ->G ->T]
%
% Nn can also be an M struct with mut and cov fields (from load_all_mutation_data2.m)
%
% C is list of best category sets discovered
% Sc is stats of each best set of categories
% Sv is stats of "vanilla" category set of (CpG transit, other C:G transit, C:G transver, A:T mut)
%
% based on category_discovery.m, Mike Lawrence August 2008
% packaged into current function 2009-12-11
% modified to accept M struct, 2010-12-01

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'method','best');   %  "greedy" or "best"
P = impose_default_value(P,'max_k',5);
P = impose_default_value(P,'mutcategs_report_filename',[]);

if isstruct(Nn)     % input is an M-type struct
  error('Not supported');
end

if size(Nn,1)~=32, error('input must have 32 rows'); end
if size(Nn,2)~=5, error('input must have 5 columns (N A C G T)'); end

orig_N = [repmat(Nn(1:16,1)',3,1);repmat(Nn(17:32,1)',3,1)];
orig_n = [Nn(1:16,[3 4 5])';Nn(17:32,[2 4 5])'];

if ~isempty(P.mutcategs_report_filename)
  report_fh = fopen(P.mutcategs_report_filename,'wt');
end

% reformat N and n (6x16) into N and n (4x4x2x3)

% dimensions:
%   (1)  5' base (1234=ACGT)
%   (2)  3' base (1234=ACGT)
%   (3)  "from" base (12=AC)
%   (4)  "to" outcome (1=transition, 2=flip_transversion, 3=skew_transversion)

n = zeros(4,4,2,3);
N = zeros(4,4,2,3);

for base5 = 1:4
  for base3 = 1:4
    for oldbase = 1:2
      for muttype = 1:3
        col = 4*(base5-1)+base3;
        switch oldbase
          case 1 % A
            switch muttype
              case 1, row = 2;  % A->G transition
              case 2, row = 3;  % A->T transversion(flip)
              case 3, row = 1;  % A->C transversion(skew)
            end
          case 2 % C
            switch muttype
              case 1, row = 6;  % C->T transition
              case 2, row = 5;  % C->G transversion(flip)
              case 3, row = 4;  % C->A transversion(skew)
            end
        end
        n(base5,base3,oldbase,muttype) = orig_n(row,col);
        N(base5,base3,oldbase,muttype) = orig_N(row,col);
end,end,end,end

Ntot = sum(N(:));
unsplit = {[1:4],[1:4],[1:2],[1:3]};

C = cell(P.max_k,1);   % output of category sets
Sc = cell(P.max_k,1);    % output of stats about category sets

% test "vanilla" set for comparison
if nargout>=3
  vanilla_leaf = [4687 4799 25343 29183];  % CpG transitions, other C:G transitions, C:G transversions, A:T mutations
  Sv = [];
  Sv.H = entropy_by_parts2(vanilla_leaf);
  subfprintf('VANILLA CATEGORY SET:\n');
  [tmp stats] = reportrule2(vanilla_leaf);
end

if ~strcmpi(P.method,'best'), fprintf('WARNING: methods other than "best" do not properly route output into C and H'); end

switch(lower(P.method))

  case 'greedy'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % navigate breakdown tree (greedy algorithm)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  initial = {unsplit};
  subfprintf('k=1:\n\n');
  H_initial = entropy_by_parts(initial);
  stats = reportrule(initial);
  
  current = initial;
  for sz=2:P.max_k
    subfprintf('k=%d:  ', sz);
    H_current = entropy_by_parts(current);
    best_new = {};
    best_dH = 0;
    for t=1:length(current)
      rest = current(setdiff(1:length(current),t));
      parent = current{t};
      for d=1:4
        tobreak = parent{d};
        p = powerset(tobreak);
        minel = min(tobreak);
        for i=2:length(p)-1
          if ismember(minel,p{i})
            child1 = parent;
            child2 = parent;
            child1{d} = p{i};
            child2{d} = setdiff(tobreak, p{i});
            new_parts = [rest;{child1};{child2}];
            H_new = entropy_by_parts(new_parts);
            dH = H_new - H_current;
            if dH<best_dH
              best_dH = dH;
              best_new = new_parts;
    end,end,end,end,end
     
    stats = reportrule(best_new);
    current = best_new;
  end % next sz

    
  case 'best-old1'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find best possible category set of a given size
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  initial = {unsplit};
  leaves = {initial};

  for sz=1:P.max_k
    subfprintf('k=%d:  ', sz);
    
    if sz>1
      old_leaves = leaves;
      leaves = {};
      for l = 1:length(old_leaves)
        leaf = old_leaves{l};
        for t=1:length(leaf)
          rest = leaf(setdiff(1:length(leaf),t));
          parent = leaf{t};
          for d=1:4
            tobreak = parent{d};
            p = powerset(tobreak);
            minel = min(tobreak);
            for i=2:length(p)-1
              if ismember(minel,p{i})
                child1 = parent;
                child2 = parent;
                child1{d} = p{i};
                child2{d} = setdiff(tobreak, p{i});
                new_leaf = [rest;{child1};{child2}];
                leaves = [leaves;{new_leaf}];
      end,end,end,end,end
      leaves = remove_duplicate_leaves(leaves);
    end
    
    best_l = 0;
    best_H = Inf;
    for l=1:length(leaves)
      H = entropy_by_parts(leaves{l});
      if H<best_H
        best_l = l;
        best_H = H;
      end
    end
    
    % report results
    best_leaf = leaves{best_l};
    stats = reportrule(best_leaf);
  end % next sz


  case 'best-old2'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find best possible category set of a given size
  % OPTIMIZED VERSION
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   initial = {unsplit};
   leaves = {initial};
   
   estimated_growth_factor = 30;   % actual split closer to 20, but let's err on generous side
   for sz=1:P.max_k
     subfprintf('k=%d:  ', sz);
     if sz>1
       old_leaves = leaves;
       typical_leaf = repmat({unsplit},sz,1);
       leaves = repmat(typical_leaf,length(old_leaves)*estimated_growth_factor,1);
       lastleaf = 0;
       for l = 1:length(old_leaves)
         leaf = old_leaves{l};
         for t=1:length(leaf)
           rest = leaf(setdiff(1:length(leaf),t));
           parent = leaf{t};
           for d=1:4
             tobreak = parent{d};
             p = powerset(tobreak);
             p = p(2:end-1);   % remove empty and full sets
             np = length(p);
             new_leaf = [rest;{parent};{parent}];
             new_leaves = repmat({new_leaf},np,1);
             child1i = length(new_leaf);
             child2i = child1i-1;
             for i=1:np
               new_leaves{i}{child1i}{d} = p{i};
               new_leaves{i}{child2i}{d} = setdiff(tobreak,p{i});
             end
             leaves(lastleaf+1:lastleaf+length(new_leaves)) = new_leaves;
             lastleaf = lastleaf + length(new_leaves);
           end,end,end
           leaves = remove_duplicate_leaves(leaves(1:lastleaf));
     end

     best_l = 0; best_H = Inf;
     for l=1:length(leaves)
       H = entropy_by_parts(leaves{l});
       if H<best_H, best_l = l; best_H = H; end
     end
     
     % report results
     best_leaf = leaves{best_l};
     stats = reportrule(best_leaf);
   end % next sz

   
  case 'best'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % superfast version based on all-integer processing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   unsplit = categ_to_int({[1 2 3 4] [1 2 3 4] [1 2] [1 2 3]});
   mask = uint16(15*16.^[0:3]);
   nybs = uint16((1:14)'*(16.^[0:3]));

   initial = unsplit;
   leaves = initial;
   for k=1:P.max_k
       subfprintf('k=%d\n',k);
       if k>1
         old_leaves = leaves;
         leaves = zeros(34^(k-1),k,'uint16');
         lastleaf = 0;
         for l=1:length(old_leaves)
           leaf = old_leaves(l,:);
           for c=1:k-1   % choose which category to split
             parent = leaf(c);
             for d=1:4   % choose which dimension to split along
               tobreak = bitand(parent,mask(d));
               tokeep = parent-tobreak;
               frags1 = bitand(tobreak,nybs(:,d));
               frags2 = tobreak-frags1;
               frags = [frags1 frags2];
               frags(~frags1|~frags2,:)=[];
               children = tokeep+frags;
               nc = size(children,1);
               leaves(lastleaf+1:lastleaf+nc,1:k-1) = repmat(leaf,nc,1);
               leaves(lastleaf+1:lastleaf+nc,[c k]) = children;
               lastleaf = lastleaf + nc;
       end,end,end
       leaves = unique(sort(leaves(1:lastleaf,:),2),'rows');
     end

     best_l = 0; best_H = Inf;

     if k<4
       for l=1:length(leaves)
         H = entropy_by_parts2(leaves(l,:));
         if H<best_H, best_l = l; best_H = H; end
       end
     else
       if k==4
         mx = (2^16)-1;
         lookup = nan(mx,1);
         for i=1:mx
           if any(bitget(i,[11 12 16])), continue; end  % those bits have no meaning
           c = int_to_categ(i);
           ns = sumpart(n,c);
           Ns = sumpart(N,c);
           f = Ns / Ntot;
           H_part = entropy(ns/Ns);
           lookup(i) = f*H_part;
         end
       end
       for l=1:length(leaves)
         H = 0;
         for i=1:k
           H = H + lookup(leaves(l,i));
         end  
         if H<best_H, best_l = l; best_H = H; end
       end
     end
           
     % report results
     best_leaf = leaves(best_l,:);
     [C{k} stats] = reportrule2(best_leaf);
     Sc{k}.H = best_H;
   end
   
 otherwise
  error('Unknown P.method');
end

if ~isempty(P.mutcategs_report_filename)
  fclose(report_fh);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions of find_mut_categs

  function subfprintf(str,varargin)
    fprintf(str,varargin{:});
    if ~isempty(P.mutcategs_report_filename)
      fprintf(report_fh,str,varargin{:});
    end
  end

  function H = entropy_by_parts(parts)
    H = 0;
    for pp=1:length(parts)
      ns = sumpart(n,parts{pp});
      Ns = sumpart(N,parts{pp});
      f = Ns / Ntot;
      H_part = entropy(ns/Ns);
      H = H + f*H_part;
    end
  end

  function H = entropy_by_parts2(partsi)
    H = 0;
    for pp=1:length(partsi)
      part = int_to_categ(partsi(pp));
      ns = sumpart(n,part);
      Ns = sumpart(N,part);
      f = Ns / Ntot;
      H_part = entropy(ns/Ns);
      H = H + f*H_part;
    end
  end

  function x = sumpart(m,cut)
    x = m(cut{1},cut{2},cut{3},cut{4});
    x = sum(x(:));
  end

  function H = entropy(p)
    if p==0 || p==1
      p1 = 0;
      p2 = 0;
    elseif p<0 || p>1
      error('p must be between zero and one');
    else
      p1 = p*log2(p);
      p2 = (1-p)*log2(1-p);
    end
    H = -(p1+p2);
  end

  function s = convert_parts_to_rule(parts)
    bases='ACGT';
    change='tfs';
    np = length(parts);
    s = cell(np,1);
    for p=1:np
      part = parts{p};
      s{p} = [bases(part{1}) '[' bases(part{3}) '->' change(part{4}) ']' bases(part{2})];
    end
  end

  function s = rulestats(parts)
    np=length(parts);
    s.n = zeros(np,1);
    s.N = zeros(np,1);
    for p=1:np
      s.N(p) = sumpart(N,parts{p});
      s.n(p) = sumpart(n,parts{p});
    end
    s.rate = s.n./s.N;
    rate_tot = fullsum(n)/fullsum(N);
    s.relrate = s.rate / rate_tot;
  end

  function stats = reportrule(parts)
    rules = convert_parts_to_rule(parts);
    rulenames = find_good_names_for_mutation_categories(rules);
    stats = rulestats(parts);
    [tmp ord] = sort(stats.relrate,'descend');
    for j=1:length(parts)
      i=ord(j);
      subfprintf('%25s   n %5.0f N %10.0f  rate %.2e (%sx)\n',...
              rulenames{i},stats.n(i),stats.N(i),stats.rate(i),...
              format_number(stats.relrate(i),3,4));
    end
    subfprintf('\n');
  end

  function [ck, stats] = reportrule2(partsi)
    for i=1:length(partsi), parts{i,1} = int_to_categ(partsi(i)); end
    rules = convert_parts_to_rule(parts);
    stats = rulestats(parts);
    tmp = parse(rules,'(\S+)\[(\S+)\->(\S+)\](\S+)',{'left','from','change','right'});
    ck = merge_structs({tmp,stats});
    ck.autoname = rules;
    ck.name = find_good_names_for_mutation_categories(ck.autoname);
    ck.type = repmat({'point'},slength(ck),1);
    [tmp ord] = sort(stats.relrate,'descend');
    for j=1:length(parts)
      i=ord(j);
      subfprintf('%25s   n %5.0f N %10.0f  rate %.2e (%sx)\n',...
              ck.name{i},stats.n(i),stats.N(i),stats.rate(i),...
              format_number(stats.relrate(i),3,4));
    end
    ck = reorder_struct(ck,ord);
    subfprintf('\n');
  end

  function LL = remove_duplicate_leaves(L)
    LL = L;
    % first convert all categories to integers
    nl = length(L);
    len = zeros(nl,1);
    for l=1:nl
      leaf = L{l};
      nc = length(leaf);
      ileaf = zeros(nc,1);
      for c=1:nc
        ileaf(c) = categ_to_int(L{l}{c});
      end
      L{l} = ileaf;
      len(l) = nc;
    end
    % compare leaves in groups of same-length leaves
    deleted = false(nl,1);
    [u ui uj] = unique(len);
    for i=1:length(u)
      idx = find(uj==i);
      x = zeros(length(idx),u(i));
      for j=1:length(idx)
        x(j,:) = L{idx(j)};
      end
      x = sort(x,2);   % sort each row
      [v vi vj] = unique(x, 'rows');
      % mark duplicates for deletion
      deleted(setdiff(idx,idx(vi))) = true;
    end
    % delete duplicates
    LL = LL(find(~deleted));
  end

  function i = categ_to_int(c)
    x = [c{1} c{2}+4 c{3}+8 c{4}+12]-1;
    i = sum(bitshift(1,x));
  end

  function c = int_to_categ(i)
    x = bitget(i,1:16);
    b = 1:4;
    c = { b(x(1:4)>0) b(x(5:8)>0) b(x(9:12)>0) b(x(13:16)>0) };
  end  

end % of find_mut_categs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ss = format_number(values,sigfigs,width)
%
% format_number(value,sigfigs,width)
%
% formats the number specified in "value" to fit in
% a character field of "width" characters, while
% displaying "sigfigs" significant figures.
%
% if decimal expansion will fit, e.g. 0.0023, then this is used;
% otherwise scientific notation is used.
%
% Mike Lawrence 2008-05-20
%
% 2011-05-17 vectorized

if sigfigs<1, error('sigfigs must be at least 1'); end

if ~isnumeric(values), error('values must be numeric'); end
if ndims(values)>2, error('N-D case not implemented'); end

ni = size(values,1); nj = size(values,2);

if ni>1e5 || nj>1e5
  fprintf('format_number: table too big, using simpler method instead\n');
  ss = cellstr(num2str(values));

else

  sigfigs = sigfigs+1;

  ss = cell(ni,nj);
  for i=1:ni, for j=1:nj
      value = values(i,j);
      
      if value>0, pdz = ceil(-log10(value))-1;    % post-decimal zeroes
      elseif value<0, pdz = ceil(-log10(-value))-1;
      else, pdz=0; % value==0
      end

      s1 = sprintf(['%.' num2str(round(sigfigs)+pdz-(abs(value)<1)) 'f'], value);
      if str2double(s1)==0 && value~=0
        s1 = sprintf(['%.' num2str(1+round(sigfigs)+pdz-(abs(value)<1)) 'f'], value);
      end

      s2 = sprintf(['%.' num2str(round(sigfigs)-2) 'd'], value);
      
      if length(s1) <= width, s = s1;
      elseif length(s2) <= length(s1), s = s2;
      else, s = s1;
      end
      
      ss{i,j} = s;
  end,end

  if ni==1 && nj==1, ss = ss{1,1}; end   % to preserve original non-vectorized behavior

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = fullsum(m)
  if length(m)==1, s = m;
  else s = fullsum(sum(m)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function names = generate_192_categ_names

bases = 'ACGT';
names = cell(192,1);
i=1;
for from=1:4
  for left=1:4
    for right=1:4
      for to=1:4
        if from==to, continue; end
        names{i,1} = [bases(left) '(' bases(from) '->' bases(to) ')' bases(right)];
        i=i+1;
end,end,end,end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function categ_list = generate_categ_context65_names()

categ_list.num = (1:65)';
x = {'A';'C';'G';'T'};
y = {}; for i=1:length(x), y = [y; regexprep(x,'^(.*)$',[x{i} '_$1'])]; end
z = {}; for i=1:length(x), z = [z; regexprep(y,'^(.*)$',[x{i} ' in $1'])]; end
categ_list.name = [z;'any N'];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function keep = grepi(pattern,strings,flag)
% case-insensitive grep

if ~exist('flag','var'), flag=0; end

keep = grep(upper(pattern),upper(strings),1);
if ~flag, keep = strings(keep); end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res,resi]=grep(reg_exp,strs,res_is_idx)
% grep

if ischar(strs)
  strs=cellstr(strs);
end

resi=find(~cellfun('isempty',cat(1,regexp(strs,reg_exp))));
res=strs(resi);

if nargout==0
  disp(catn(res));
end

if exist('res_is_idx','var') && res_is_idx
  [res,resi]=exchange_vars(res,resi);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function keep = grepv(pattern,strings,flag)
% inverse grep

if ~exist('flag','var'), flag=0; end

toss = grep(pattern,strings,1);
keep = setdiff((1:length(strings))',toss);
if ~flag, keep = strings(keep); end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = hist2d(y,x,ybins,xbins)

y = as_column(y);
x = as_column(x);
if length(y)~=length(x), error('y and x need to be same length'); end
h = zeros(length(ybins),length(xbins));
for i=1:length(y)
  yidx = find(y(i)>=ybins,1,'last');
  if isempty(yidx), yidx=length(ybins); end
  xidx = find(x(i)>=xbins,1,'last');
  if isempty(xidx), xidx=length(xbins); end
  h(yidx,xidx)=h(yidx,xidx)+1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S2 = keep_fields(S,flds)
% keep_fields(S,flds)
%
% given struct <S> and cell-array-of-strings <flds>,
% returns struct <S2> which has only those fields specified.
%
% Mike Lawrence 2009-04-24

if ischar(flds), flds = {flds}; end

S2=[];
for i=1:length(flds)
  if isempty(S)
    f = [];
  else
    f = getfield(S,flds{i});
  end
  S2=setfield(S2,flds{i},f);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = list2cell(varargin)

for i=1:nargin
  c{i,1} = varargin{i};
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = mapacross(A,B,C,filler)
% D = mapacross(A,B,C[,filler])
%
% for each item in A,
%   looks for that item in B.
%   if found, returns the corresponding item in C;
%   if not found, returns "filler".
% returns the items as D.
%
% Mike Lawrence 2009-07-14

if nargin~=3 && nargin~=4
  error('usage: D = mapacross(A,B,C[,filler])');
end

if ~isnumeric(A) && ~iscell(A) && ~islogical(A)
  error('A should be numeric, cell, or logical');
end
if ~isnumeric(B) && ~iscell(B) && ~islogical(B)
  error('B should be numeric, cell, or logical');
end
if ~isnumeric(C) && ~iscell(C) && ~islogical(C)
  error('C should be numeric, cell, or logical');
end

if length(B) ~= length(C), error('length(B) ~= length(C)'); end

idx = listmap(A,B);
idx2 = find(~isnan(idx));

if ndims(C)>2 || (size(C,1)>1 && size(C,2)>1)
  flag=true;
  dsz = size(C); dsz(1) = length(A);
else
  flag=false;
  dsz = [length(A) 1];
end

if iscell(C)
  if ~exist('filler','var'), filler = ''; end
  D = repmat({filler},dsz);
elseif isnumeric(C)
  if ~exist('filler','var'), filler = nan; end
  D = filler*ones(dsz);
elseif islogical(C)
  if length(idx2)~=length(idx), error('Missing values when using logical output: behavior undefined'); end
  D = false(dsz);
else
  error('unknown output type');
end

if flag
  if ndims(C)>7, error('>7-D not supported'); end
  D(idx2,:,:,:,:,:,:) = C(idx(idx2),:,:,:,:,:,:);
else
  D(idx2) = C(idx(idx2));
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = merge_structs(slist)
%
% merge_structs(slist)
%
% slist should be a cell array of structs,
%   all of which have the slength
% 
% returns a new struct which has all the fields of the input structs
%
% Mike Lawrence 2009-02-04
%

if ~iscell(slist)
  error('input should be a cell array of structs');
end

S=[];
for i=1:numel(slist)
  if ~isstruct(slist{i})
    error('input should be a cell array of structs');
  end
  s = slist{i};
  if i==1
    S = slist{i};
  else
    if slength(s) ~= slength(S)
      error('all input structs must have same length');
    end
    fields = fieldnames(s);
    for fno=1:length(fields), fn = fields{fno};
      f = getfield(s, fn);
%      if isfield(S,fn), fprintf('Warning: duplicate field %s overwritten\n',fn);end
      S = setfield(S,fn,f);
    end
  end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u ui uj] = nanunique(a,varargin)

[u ui uj] = unique(a,varargin{:});

if isnumeric(a)
  idx = find(isnan(u));
  if length(idx)>1
    if idx(end)~=length(u) || any(diff(idx)~=1), error('expected unique to place all NaNs at end'); end
    uj(isnan(a)) = idx(1);
    u(idx(2:end)) = [];
    ui(idx(2:end)) = [];
  end
end
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = parse_in(x,parsefield,pattern,newfields,numidx)

if ~exist('numidx','var'), numidx=[]; end

tmp = parse(getfield(x,parsefield),pattern,newfields,numidx);
x = merge_structs({x,tmp});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = powerset(S)
  L = length(S);
  P = cell(2^L,1);
  for i=0:(2^L)-1
    idx = find(bitand(i,2.^(0:L-1)));
    P{i+1}=S(idx);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = rdir(rootdir,varargin)
% Lists the files in a directory and its sub directories. 
% 
% function [D] = rdir(ROOT,TEST)
%
% Recursive directory listing.
%
% ROOT is the directory starting point and includes the 
% wildcard specification.
% The function returns a structure D similar to the one 
% returned by the built-in dir command. 
% There is one exception, the name field will include 
% the relative path as well as the name to the file that 
% was found.
% Pathnames and wildcards may be used. Wild cards can exist
% in the pathname too. A special case is the double * that
% will match multiple directory levels, e.g. c:\**\*.m. 
% Otherwise a single * will only match one directory level.
% e.g. C:\Program Files\Windows *\
%
% TEST is an optional test that can be performed on the 
% files. Two variables are supported, datenum & bytes.
% Tests are strings similar to what one would use in a 
% if statement. e.g. 'bytes>1024 & datenum>now-7'
%
% If not output variables are specified then the output is 
% sent to the screen.
%
% See also DIR
%
% examples:
%   D = rdir('*.m');
%     for ii=1:length(D), disp(D(ii).name); end;
%
%   % to find all files in the current directory and sub directories
%   D = rdir('**\*')
%
%   % If no output is specified then the files are sent to 
%   % the screen.
%   rdir('c:\program files\windows *\*.exe');
%   rdir('c:\program files\windows *\**\*.dll');
%
%   % Using the test function to find files modified today
%   rdir('c:\win*\*','datenum>floor(now)');
%   % Using the test function to find files of a certain size
%   rdir('c:\program files\win*\*.exe','bytes>1024 & bytes<1048576');
%

% use the current directory if nothing is specified
if ~exist('rootdir','var'),
  rootdir = '*';
end;

% split the file path around the wild card specifiers
prepath = '';       % the path before the wild card
wildpath = '';      % the path wild card
postpath = rootdir; % the path after the wild card
I = find(rootdir==filesep,1,'last');
if ~isempty(I),
  prepath = rootdir(1:I);
  postpath = rootdir(I+1:end);
  I = find(prepath=='*',1,'first');
  if ~isempty(I),
    postpath = [prepath(I:end) postpath];
    prepath = prepath(1:I-1);
    I = find(prepath==filesep,1,'last');
    if ~isempty(I),
      wildpath = prepath(I+1:end);
      prepath = prepath(1:I);
    end;
    I = find(postpath==filesep,1,'first');
    if ~isempty(I),
      wildpath = [wildpath postpath(1:I-1)];
      postpath = postpath(I:end);
    end;
  end;
end;
% disp([' "' prepath '" ~ "' wildpath '" ~ "' postpath '" ']);

if isempty(wildpath),
  % if no directory wildcards then just get file list
  D = dir([prepath postpath]);
  keep = true(length(D),1);
  for ii = 1:length(D)
    if strcmp(D(ii).name,'.'), keep(ii) = false; end
    if strcmp(D(ii).name,'..'), keep(ii) = false; end
  end
  D = D(keep);
%  D([D.isdir]==1) = [];
  for ii = 1:length(D),
%    if (~D(ii).isdir),
      D(ii).name = [prepath D(ii).name];
%    end;
  end;

  % disp(sprintf('Scanning "%s"   %g files found',[prepath postpath],length(D)));

elseif strcmp(wildpath,'**'), % a double wild directory means recurs down into sub directories

  % first look for files in the current directory (remove extra filesep)
  D = rdir([prepath postpath(2:end)]);

  % then look for sub directories
  Dt = dir(''); 
  tmp = dir([prepath '*']);
  % process each directory
  for ii = 1:length(tmp),
    if (tmp(ii).isdir && ~strcmpi(tmp(ii).name,'.') && ~strcmpi(tmp(ii).name,'..') ),
      Dt = [Dt; rdir([prepath tmp(ii).name filesep wildpath postpath])];
    end;
  end;
  D = [D; Dt];

else
  % Process directory wild card looking for sub directories that match
  tmp = dir([prepath wildpath]);
  D = dir(''); 
  % process each directory found
  for ii = 1:length(tmp),
    if (tmp(ii).isdir && ~strcmpi(tmp(ii).name,'.') && ~strcmpi(tmp(ii).name,'..') ),
      D = [D; rdir([prepath tmp(ii).name postpath])];
    end;
  end;
end;


% Apply filter
if (nargin>=2 && ~isempty(varargin{1})),
  date = [D.date];
  datenum = [D.datenum];
  bytes = [D.bytes];

  try
    eval(sprintf('D((%s)==0) = [];',varargin{1})); 
  catch
    warning('Error: Invalid TEST "%s"',varargin{1});
  end;
end;

% display listing if no output variables are specified
if nargout==0,
  pp = {'' 'k' 'M' 'G' 'T'};
  for ii=1:length(D), 
    sz = D(ii).bytes;
    if sz<=0,
      disp(sprintf(' %31s %-64s','',D(ii).name)); 
    else
      ss = min(4,floor(log2(sz)/10));
      disp(sprintf('%4.0f %1sb   %20s   %-64s ',sz/1024^ss,pp{ss+1},D(ii).date,D(ii).name)); 
    end;
  end;
else
  % send list out
  varargout{1} = D;
end;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function q = remove_ansi(q)
  q = regexprep(q,[char(27) '\[[^m]*m'],'');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function require_fields(T,fields)

  if ~iscell(fields)
    fields = {fields};
  end

  for i=1:length(fields)
    if ~isfield(T,fields{i})
      error(['Structure is missing required field "' fields{i} '"']);
    end
  end
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = rmfield_if_exist(S,f)
  if ~iscell(f), f = {f}; end
  for i=1:length(f)
    if isfield(S,f{i}), S = rmfield(S,f{i}); end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_textfile(t,filename)
  out = fopen(filename,'wt');
  fwrite(out,t,'char');
  fclose(out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = stringsplice(A,dim,sep)
% stringsplice(A,dim,sep)
%
% given a cell matrix of strings <A>
% and a dimension <dim> (default=1),
% concatenates strings along the dimension and returns a cell array of strings.
%
% Mike Lawrence 2009-03-03

if nargin==2
  if ischar(dim)
    sep=dim;
    clear dim;
  end
end

if nargin==3
  if ischar(dim) && isnumeric(sep)
    tmp = dim;
    dim = sep;
    sep = tmp;
  elseif ischar(dim) || isnumeric(sep)
    error('unable to parse input parameters');
  end
end

if ~exist('dim','var'), dim=1; end
if ~exist('sep','var'), sep=''; end

if dim==1, B = cell(size(A,1),1);
else B = cell(1,size(A,2)); end

for i=1:length(B), if ~mod(i,1e5), fprintf('%d/%d ',i,length(B)); end
  if dim==1, B{i} = concat(A(i,:),sep);
  else B{i} = concat(A(:,i),sep); end
end, if i>=1e5, fprintf('\n'); end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ct,ua,ub] = xcount(a, b, sort_by_frequency, truncate_length)
% xcount(a,b,sort_by_frequency)
%
% Finds the unique pairs across two arrays (a and b) and counts how many of each pair there
%
% If sort_by_frequency is 1 or -1 (reverse), then table is sorted.
%
% if truncate_length is specified, trims strings to specified length (to improve visual fit)
%
% Mike Lawrence 2008-06-16
% 2012-12-12 speedup for working with numbers

if size(a,1)==1, a=a'; end
if size(b,1)==1, b=b'; end

if length(a)~=length(b), error('Arrays must be same length'); end

% find unique elements and make table

[ua uai uaj] = nanunique(a);
nua = length(ua);

[ub ubi ubj] = nanunique(b);
nub = length(ub);

% make sure ua and ub are cell arrays of strings
% (if they aren't, convert them)

if ~iscell(ua)
  tmp=ua;
  ua = cell(length(tmp),1);
  for i=1:length(tmp)
    if ischar(tmp(i))
      ua(i) = {tmp(i)};
    else
      ua(i) = {num2str(tmp(i))};
end,end,end

if ~iscell(ub)
  tmp=ub;
  ub = cell(length(tmp),1);
  for i=1:length(tmp)
    if ischar(tmp(i))
      ub(i) = {tmp(i)};
    else
      ub(i) = {num2str(tmp(i))};
end,end,end

% if truncate_length specified, trim values
if exist('truncate_length','var')
  for i=1:length(ua)
    ua{i} = ua{i}(1:min(truncate_length*2,length(ua{i})));
  end
  for i=1:length(ub)
    ub{i} = ub{i}(1:min(truncate_length,length(ub{i})));
  end
end

% compute 2D histogram
ct = zeros(nua,nub);
for aj=1:nua
  for bj=1:nub
    ct(aj,bj) = sum(uaj==aj & ubj==bj);
  end
end

% compute totals
tota = sum(ct,2);
totb = sum(ct,1);

% sort if requested (rows by tota, cols by totb)
if exist('sort_by_frequency', 'var') && sort_by_frequency
  if sort_by_frequency == 1
    [tmp orda] = sort(tota);
    [tmp ordb] = sort(totb);
  elseif sort_by_frequency == -1
    [tmp orda] = sort(tota, 'descend');
    [tmp ordb] = sort(totb, 'descend');
  end
  ct = ct(orda,ordb);
  tota = tota(orda);
  totb = totb(ordb);
  ua = ua(orda);
  ub = ub(ordb);
end

% make cell table, including totals

tbl = cell(nua+2,nub+2);
for aj=1:nua
  for bj=1:nub
    tbl{aj+1,bj+1} = num2str(ct(aj,bj));
  end
end

for aj=1:nua
  tbl{aj+1,1} = ua{aj};
  tbl{aj+1,end} = num2str(tota(aj));
end

for bj=1:nub
  tbl{1,bj+1} = ub{bj};
  tbl{end,bj+1} = num2str(totb(bj));
end

tbl{end,1} = 'TOTAL';
tbl{1,end} = 'TOTAL';
tbl{1,1} = '';
tbl{end,end} = num2str(sum(tota));

% measure columns

L=zeros(nua+2,nub+2);
for aj=1:nua+2
  for bj=1:nub+2
    L(aj,bj)=length(tbl{aj,bj});
  end
end

W=zeros(nub+2,1);
for bj=1:nub+2
  W(bj) = max(L(:,bj));
end

if nargout==0
  % print table
  fprintf('\n');
  for aj=1:nua+2
    if aj==nua+2, fprintf('\n'); end
    fprintf('    ');
    for bj=1:nub+2
      if bj==2, fprintf('  '); end
      if bj==nub+2, fprintf('    '); end
      fprintf(['%' num2str(W(bj)) 's  '], tbl{aj,bj});
    end
    fprintf('\n');
  end
  fprintf('\n');
end

if nargout<3, clear ub; end
if nargout<2, clear ua; end
if nargout<1, clear ct; end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s=reorder_struct_exclude(s,order)
  if islogical(order), order = find(order); end
  s = reorder_struct(s,setdiff(1:slength(s),order));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



