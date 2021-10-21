//for the question 6

exec('Genomics_signal_proc1.sce');
exec('Genomics_signal_proc2.sce');


//accession list
file_list = ['NZ_CP066047.1','NZ_CP028172.1','NZ_CP030194.1','NZ_CP030231.1','NZ_CP030238.1','NZ_CP037891.1','NZ_CP040380.1','NZ_CP046277.1','NZ_CP046279.1','NZ_CP046280.1','NZ_CP046291.1','NZ_CP053581.1','NZ_CP060508.1','NZ_CP069518.1','NZ_CP075108.1']
list_size = size(file_list,2);

//Iterate over all the files
for i=1:list_size
    disp(file_list(i))
    [gp,gn,ncp,ncn]= get_protein_pos_array(file_list(i)+'.txt');//('NZ_CP060508.1.txt');//
    
    fasta_in =  file_list(i)+'.fasta ';//'NZ_CP060508.1.fasta';//
    
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


    //50 base check
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


    test_set_sense = methionine_checked_sense;
    test_set_antisense = methionine_checked_antisense;

    mprintf('Iteration for threshold starts here \n ------------------------- \n', consensus_score_reduced);

    //mprintf('Total valid sense genes in test set  = %i \n', sum(test_set_sense));
    //mprintf('Total valid antisense genes in test set  = %i \n', sum(test_set_antisense));

    num_of_sense_genes = sum(test_set_sense);
    num_of_antisense_genes = sum(test_set_antisense);

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
        /*
        mprintf('Total valid sense genes in test set  = %i \n', num_of_sense_genes);

        mprintf('Total detected promoters in sense using original ppm  = %i | .. %f  |threshold = %i \n', detected_promoters, detected_promoters/num_of_sense_genes*100 , consensus_thresh);

        mprintf('Total detected promoters in sense using reduced ppm  = %i | .. %f  | threshold = %i \n', detected_promoters_reduced, detected_promoters_reduced/num_of_sense_genes*100, consensus_thresh);

        mprintf('---\n')
        */
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

        /*
        mprintf('Total valid antisense genes in test set  = %i \n', num_of_antisense_genes);

        mprintf('Total detected promoters in anti-sense using original ppm  = %i | .. %f  | threshold = %i \n', detected_promoters_antisense,detected_promoters_antisense/num_of_antisense_genes*100, consensus_thresh);

        mprintf('Total detected promoters in anti-sense using reduced ppm  = %i | .. %f  | threshold = %i \n', detected_promoters_reduced_antisense,detected_promoters_reduced_antisense/num_of_antisense_genes*100, consensus_thresh);

        mprintf('---\n')
        */
        total_genes = num_of_antisense_genes+num_of_sense_genes;
        detected = detected_promoters_antisense+detected_promoters;
        detected_reduced = detected_promoters_reduced_antisense +detected_promoters_reduced;

        mprintf('--Total valid genes in test set  = %i \n', total_genes);

        mprintf('--Total detected promoters using original ppm  = %i | .. %f  | threshold = %i \n', detected,detected/total_genes*100, consensus_thresh);

        mprintf('--Total detected promoters in using reduced ppm  = %i | .. %f  | threshold = %i \n', detected_reduced,detected_reduced/total_genes*100, consensus_thresh);


        mprintf('-------------------------------------------------------------\n')
    end

end
