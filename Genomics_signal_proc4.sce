//Question 5 part 2

entropy_thresh = 0.05;
consensus_thresh = -1;
mprintf('Used entropy threshold = %f \n', entropy_thresh);

// Get consensus score 
//use the correct consensus sequence here after checking with part 3 code
seq_con = ascii('ATTATAAAAA')//ascii('ATTAAAACAG' ); 

consensus_score=stat_align(seq_con,pribnow_ppm);
mprintf('Consensus score = %f \n', consensus_score);

consensus_score_reduced=max(stat_align_entropy(seq_con,pribnow_ppm,entropy_thresh));
mprintf('Reduced Consensus score = %f \n', consensus_score_reduced);

mprintf('Iteration for threshold starts here \n ------------------------- \n', consensus_score_reduced);

//mprintf('Total valid sense genes in test set  = %i \n', sum(test_set_sense));
//mprintf('Total valid antisense genes in test set  = %i \n', sum(test_set_antisense));

for consensus_thresh =-1:-1:-5
    detected_promoters = 0;
    detected_promoters_reduced = 0;

    detected_promoters_antisense = 0;
    detected_promoters_reduced_antisense = 0;
    for n = 1:sense_genes_N
        if (test_set_sense(n))
            sequence = get_fasta_at(fasta_in, gp(n,1)-pribnow_start,gp(n,1)-pribnow_stop,1);
            
            score = stat_align(sequence,pribnow_ppm);
            score_reduced = stat_align_entropy(sequence,pribnow_ppm,entropy_thresh);

            if((max(score)-consensus_score)>consensus_thresh) then
                detected_promoters = detected_promoters +1;
            end

            if((max(score_reduced)-consensus_score_reduced)>consensus_thresh) then
                detected_promoters_reduced = detected_promoters_reduced +1;
            end
        end
    end

    mprintf('Total valid sense genes in test set  = %i \n', sum(test_set_sense));

    mprintf('Total detected promoters in sense using original ppm  = %i | threshold = %i \n', detected_promoters, consensus_thresh);

    mprintf('Total detected promoters in sense using reduced ppm  = %i | threshold = %i \n', detected_promoters_reduced, consensus_thresh);

    mprintf('---\n')

    for n = 1:antisense_genes_N
        if (test_set_antisense(n))
            sequence = get_fasta_at(fasta_in, gn(n,2)+pribnow_stop, gn(n,2)+ pribnow_start,0);
            
            score = stat_align(sequence,pribnow_ppm);
            score_reduced = stat_align_entropy(sequence,pribnow_ppm,entropy_thresh);

            if((max(score)-consensus_score)>consensus_thresh) then
                detected_promoters_antisense = detected_promoters_antisense +1;
            end

            if((max(score_reduced)-consensus_score_reduced)>consensus_thresh) then
                detected_promoters_reduced_antisense = detected_promoters_reduced_antisense +1;
            end
        end
    end

    mprintf('Total valid antisense genes in test set  = %i \n', sum(test_set_antisense));

    mprintf('Total detected promoters in anti-sense using original ppm  = %i | threshold = %i \n', detected_promoters_antisense, consensus_thresh);

    mprintf('Total detected promoters in anti-sense using reduced ppm  = %i | threshold = %i \n', detected_promoters_reduced_antisense, consensus_thresh);

    mprintf('---\n')

    total_genes = num_of_antisense_genes+num_of_sense_genes;
    detected = detected_promoters_antisense+detected_promoters;
    detected_reduced = detected_promoters_reduced_antisense +detected_promoters_reduced;

    mprintf('--Total valid genes in test set  = %i \n', total_genes);

    mprintf('--Total detected promoters using original ppm  = %i | .. %f  | threshold = %i \n', detected,detected/total_genes*100, consensus_thresh);

    mprintf('--Total detected promoters in using reduced ppm  = %i | .. %f  | threshold = %i \n', detected_reduced,detected_reduced/total_genes*100, consensus_thresh);


    mprintf('-------------------------------------------------------------\n')
end
