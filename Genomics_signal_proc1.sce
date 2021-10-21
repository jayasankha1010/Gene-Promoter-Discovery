

//Run the standalone functions 
exec('bm4321_sequence_alignment_func.sce');
exec('bm4321_gene_prom_region.sce');

//read the protien table
[gp,gn,ncp,ncn]= get_protein_pos_array('NZ_CP053581.1.txt');

//fasta file name
fasta_in =  'NZ_CP053581.1.fasta ';

sense_genes_N = size(gp,1);
antisense_genes_N = size(gn,1)

pribnow_start = 30; // Start position of pribnow Box
pribnow_stop = 5;   // Stop position of pribnow Box

disp('Number of genes in sense strand')
disp(sense_genes_N)
//disp(gp)

disp('Number of genes in anti-sense strand')
disp(antisense_genes_N)
//disp(gn)

//Question 1
//Remove genes that at the front due to miscalculations

sense_consider = ones(1:sense_genes_N);
antisense_consider = ones(1:antisense_genes_N);

for n_key=1:sense_genes_N
    if ((gp(n_key,1)-pribnow_start)<=0) then
        sense_consider(n_key) = 0;
    end
end

for n_key=1:sense_genes_N-1
    if ((gp(n_key+1,1)-gp(n_key,2))<50) then
        sense_consider(n_key+1) = 0;
    end
end


// Neglect the genes that are less than 50 bases upstream
considerable_sense_genes_N = sum(sense_consider);
disp('Number of Considerable Genes in sense strand(considering gap of 50)')
disp(considerable_sense_genes_N)

for n_key=1:antisense_genes_N-1
    if ((gn(n_key+1,1)-gn(n_key,2))<50) then
        antisense_consider(n_key) = 0;
    end
end

for n_key=1:sense_genes_N
    if ((gp(n_key,1)-pribnow_start)<=0) then
        sense_consider(n_key) = 0;
    end
end

considerable_antisense_genes_N = sum(antisense_consider);
disp('Number of Considerable Genes in anti-sense strand(considering gap of 50)')
disp(considerable_antisense_genes_N)



//Methionine check
//1. Sense strand
methionine_checked_sense = sense_consider;
for n_key=1:sense_genes_N
    m_check = get_fasta_at(fasta_in,gp(n_key,1),gp(n_key,1)+3,1);
    if(sense_consider(n_key)==1) then
        if (m_check==ascii('ATG') || m_check==ascii('GTG')) then
            //Check if first codon is Methionine
            //SvalidGenes = SvalidGenes+1;
        else
            methionine_checked_sense(n_key)=0;
        end
    end
end

number_of_valid_genes_sense = sum(methionine_checked_sense);
disp('Number of Valid Genes in sense strand(Considering ATG)')
disp(number_of_valid_genes_sense);

//2. Antisense strand
methionine_checked_antisense = antisense_consider;
for n_key=1:antisense_genes_N
    m_check = get_fasta_at(fasta_in,gn(n_key,2)-2,gn(n_key,2)+1,0);
    //disp(m_check)
    if(antisense_consider(n_key)==1) then
        if (m_check==ascii('ATG') || m_check==ascii('GTG')) then
            //Check if first codon is Methionine
            //ASvalidGenes = ASvalidGenes+1;
        else
            methionine_checked_antisense(n_key)=0;
        end
    end
end

number_of_valid_genes_antisense = sum(methionine_checked_antisense);
disp('Number of Valid Genes in antisense strand(Considering ATG)')
disp(number_of_valid_genes_antisense);

//Question 2
//Locate pribnow box, Using the first 1000 sequences obtain a position probability matrix PPM with 10 positions
// Get region with Pribnow Box


pribnow_len = pribnow_start-pribnow_stop; // Length of region with pribnow Box
pribnow_mat_len = 10; // Length of PPM of pribnow Box

pribnow_ppm_limit = 1000; //limit PPM to 1000 sequences
pribnow_skipped_pos = [];// to store sequences not included in ppm

k_init = 0.01; // Initial probability (to allow log calculations)

pribnow_pos_freq_matrix = k_init*ones(4,pribnow_mat_len); // Create empty frequency matrix for nucleotides to subsequently generate PPM of pribnow Box

pribnow_query = ascii('TATAAT'); // pribnow box of Salmonella Enterica

test_set_sense = methionine_checked_sense;
test_set_antisense = methionine_checked_antisense;

// Start search    
for n_key=1:sense_genes_N
    
    if(methionine_checked_sense(n_key)==1) then
        
        // Get position of pribnow Box
        pribnow_seq = get_fasta_at(fasta_in,gp(n_key,1)-pribnow_start,gp(n_key,1)-pribnow_stop,1); // Potential region to search
        [ax,ay,pribnow_pos] = traceback_prom(pribnow_seq,pribnow_query,1,-1,gap_penalty); // Promoter alignment match A or T (W) with A or T (W)

        if ((pribnow_pos>0)&(pribnow_pos<(pribnow_len-pribnow_mat_len))&(sum(pribnow_pos_freq_matrix,1)(1)<pribnow_ppm_limit)) then
            // Update the pribnow position frequency matrix
            pribnow_post_prom_seq = pribnow_seq((pribnow_pos+1):length(pribnow_seq));
            pribnow_pos_freq_matrix = update_pos_freq_mat(pribnow_post_prom_seq,pribnow_pos_freq_matrix,pribnow_mat_len);

            //remove this from test set
            test_set_sense(n_key) = 0;
        end
        
    end
end
disp('Test set size sense strand')
disp(sum(test_set_sense))


disp('Pribnow Position Frequency Matrix')
disp(pribnow_pos_freq_matrix)

pribnow_ppm = get_ppm(pribnow_pos_freq_matrix); // Generate pribnow Box PPM
disp('Pribnow Box PPM');
disp(pribnow_ppm); 


//Question 3
//Using suitable entrophy measure, eliminate the reedundant positions of (2). Plot the distribution of the 
//entrophy vs, number of positions and hence, select a suitable threshold 

[pribnow_w,pribnow_su]=ppm_info(pribnow_ppm,[0.25 0.25 0.25 0.25]); // Find Pribnow Box entropy
entropy = sum(pribnow_su,1);
disp(entropy);

figure;
plot(entropy);
xlabel('Position');
ylabel('Entropy');



