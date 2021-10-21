//for Q5 part 1

//Edit following value
pribnow_ppm_limit = 100; //limit PPM to 1000 sequences
mprintf('Pribnow Box PPM for limit = %i \n', pribnow_ppm_limit);

k_init = 0.01; // Initial probability (to allow log calculations)

pribnow_pos_freq_matrix = k_init*ones(4,pribnow_mat_len); // Create empty frequency matrix for nucleotides to subsequently generate PPM of pribnow Box

count=0;
// Start search    
for n=1:sense_genes_N
    n_key = floor((rand(1))*sense_genes_N)+1; //Randomly selecting genes
    if(methionine_checked_sense(n_key)==1) then
        
        // Get position of pribnow Box
        pribnow_seq = get_fasta_at(fasta_in,gp(n_key,1)-pribnow_start,gp(n_key,1)-pribnow_stop,1); // Potential region to search
        [ax,ay,pribnow_pos] = traceback_prom(pribnow_seq,pribnow_query,1,-1,gap_penalty); // Promoter alignment match A or T (W) with A or T (W)

        if ((pribnow_pos>0)&(pribnow_pos<(pribnow_len-pribnow_mat_len))&(sum(pribnow_pos_freq_matrix,1)(1)<pribnow_ppm_limit)) then
            // Update the pribnow position frequency matrix
            pribnow_post_prom_seq = pribnow_seq((pribnow_pos+1):length(pribnow_seq));
            pribnow_pos_freq_matrix = update_pos_freq_mat(pribnow_post_prom_seq,pribnow_pos_freq_matrix,pribnow_mat_len);
            count=count+1;
    
        end
    end
end
disp('count')
disp(count)

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
ylabel('Entropy')

