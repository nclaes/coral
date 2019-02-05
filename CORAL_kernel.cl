
//  USING BANDED MYERS BIT VECTOR ALGORITHM
#define READ_LENGTH RLEN 
#define W WORD_LENGTH   //W is the WORD_LENGTH
#define NUM_OF_QGRAMS QGRAMS_IN_READ
#define GENOME_LENGTH GENOME_LEN
#define ERROR PERMISSIBLE_ERROR
#define CAND_LOC_PER_READ CANDIDATES_PER_READ
#define MASK UINT_WITH_MSB_ONE // 32-bit unsigned integer with MSB as '1' rest all bits '0'
#define MIN_QG_LEN	MIN_QGRAM_LEN
#define EXTRA_BASES EXCESS_BASE
#define NO_OF_ALPHABETS 4
#define NO_OF_ALPHABETS_INC_N 5


__kernel void coral(__global uchar* genome, __global char* Reads, __global uint* SA, __global uint* tally, __global uint* F,__global uint* cand_locs_per_read, __global uchar* genomic_strand_and_ED_for_mapped_reads, __global uint* endpos_for_mapped_reads)
{	
	char RF[READ_LENGTH], RR[READ_LENGTH];  // RF - Read forward and RR - Read reverse
	uint local_F[NO_OF_ALPHABETS_INC_N];   // This should be 5 (no of alphabets + 1)
	int i=0,j=0,k=0,x=0,s=0,c = ERROR + 1;
	int band_len = READ_LENGTH + 2*ERROR, constant1 = W - 1 + READ_LENGTH - c;
	uint no_of_locations=0, B_F[NO_OF_ALPHABETS_INC_N] = {0,0,0,0,0}, B_R[NO_OF_ALPHABETS_INC_N] = {0,0,0,0,0}, B[NO_OF_ALPHABETS_INC_N] = {0,0,0,0,0};
	uint X, D0, HP, HN, VP, VN, sa_start_pos, occurences;
	int edit_dist, temp_location, verif_start_pos_in_genome, score;
	uint gid = get_global_id(0)*READ_LENGTH, last_successful_location=0;
	uchar unused_eb=0;
	ushort k_mer_using_eb=0;
	// short verif_start_pos[] = {89,73,56,40,22,5}; // (q-gram_start_position + ERROR) This is to reach the start position of the reference genome where verification will being, keeping in mind (n+2e)
	int eb=EXTRA_BASES, count = 0;
	for(i=0; i < READ_LENGTH; i++)
	{				
		switch(Reads[gid + i])
		{
			case 'T':
				RF[i] = 3;
				RR[READ_LENGTH-1-i] = 0;
				break;
			case 'C':
				RF[i] = 1;
				RR[READ_LENGTH-1-i] = 2;
				break;
			case 'G':
				RF[i] = 2;
				RR[READ_LENGTH-1-i] = 1;
				break;
			case 'A':
				RF[i] = 0;
				RR[READ_LENGTH-1-i] = 3;
				break;
			default:
				RF[i] = 4;
				RR[READ_LENGTH-1-i] = 4;  // when the base is 'N'
				break;
		}		
	}	
	for(j=0; j < c; j++) //preprocessing for banded Myers bit-vector algorithm
	{
		B_F[RF[j]] = B_F[RF[j]] | (1 << (W - c + j));
		B_R[RR[j]] = B_R[RR[j]] | (1 << (W - c + j));
	}
	
	for(i=0; i < NO_OF_ALPHABETS + 1; i++)
	{
		local_F[i] = F[i];
	}
	//Filtration testing
	if(eb > 0)
	{
		x = READ_LENGTH;	eb = EXTRA_BASES;
		for (i = 0; i < NUM_OF_QGRAMS; i++)
		{
			x--;
			occurences = 0;
			if(RF[x] != 4)
			{
				sa_start_pos = local_F[RF[x]];
				occurences = local_F[RF[x] + 1] - local_F[RF[x]];
				k = 1; x--;count = 1;			
				while(k)
				{	
					if(RF[x] == 4)			
					{
						if(count < MIN_QG_LEN)
						{
							x = x + count - MIN_QG_LEN;
						}
						occurences = 0;			
						break;					
					}
					temp_location = tally[(sa_start_pos -1)*NO_OF_ALPHABETS + RF[x]];
					occurences = tally[(sa_start_pos + occurences - 1)*NO_OF_ALPHABETS + RF[x]] - temp_location;
					sa_start_pos = local_F[RF[x]] + temp_location;
					if(count < MIN_QG_LEN)
					{
						x--;
						count++;
					}
					else if(occurences <= 1000)
					{
						k=0;
					}
					else if(occurences > 1000)
					{
						if(eb > 0 && RF[x-1] != 4)
						{
							x--;
							eb--;
							k_mer_using_eb = k_mer_using_eb | (1<<i);
						}
						else
						{
							k=0;
						}
					}
				}							
			}
			else
			{
				x = x-MIN_QG_LEN;
			}
		}
		unused_eb = eb; k_mer_using_eb = ~k_mer_using_eb;
	}
	// Filtration
	gid = get_global_id(0)*CAND_LOC_PER_READ;
	x = READ_LENGTH;	eb = EXTRA_BASES;
	for (i = 0; i < NUM_OF_QGRAMS; i++)
	{
		x--;
		occurences = 0;
		if(RF[x] != 4)
		{
			sa_start_pos = local_F[RF[x]];
			occurences = local_F[RF[x] + 1] - local_F[RF[x]];
			k = 1; x--;count = 1;			
			while(k)
			{	

				if(RF[x] == 4)			
				{
					if(count < MIN_QG_LEN)
					{
						x = x + count - MIN_QG_LEN;
					}
					occurences = 0;			
					break;					
				}
				temp_location = tally[(sa_start_pos -1)*NO_OF_ALPHABETS + RF[x]];
				occurences = tally[(sa_start_pos + occurences - 1)*NO_OF_ALPHABETS + RF[x]] - temp_location;
				sa_start_pos = local_F[RF[x]] + temp_location;
				if(count < MIN_QG_LEN)
				{
					x--;
					count++;
				}
				else if(occurences <= 1000)
				{
					if((k_mer_using_eb & (1 << i)) > 0 && (unused_eb > 0) && (eb > 0) && (RF[x-1] != 4))
					{
						unused_eb--;
						x--;
						eb--;
						k_mer_using_eb = k_mer_using_eb & (~(1 << i));
					}
					else
					{
						k=0;
					}					
				}
				else if(occurences > 1000)
				{
					if(eb > 0 && RF[x-1] != 4)
					{
						x--;
						eb--;
					}
					else
					{
						k=0;
					}
				}
			}			
			occurences = (occurences > 1000)?1000:occurences;
			// no_of_locations = no_of_locations + occurences;
			for(j = 0; j < occurences; j++)
			{	
				verif_start_pos_in_genome = SA[sa_start_pos + j] - x - ERROR; //overflow problem can occur if genome length is greater than 2^31 as this variable is signed integer. Keep in mind when dealing with whole genome		
				if(no_of_locations >= CAND_LOC_PER_READ || (verif_start_pos_in_genome > last_successful_location && verif_start_pos_in_genome < last_successful_location  + 4*ERROR))
				{
					// printf("%u 	 %u	%u 		%u	%u 		%u\n",i, x, SA[sa_start_pos + j], verif_start_pos_in_genome, verif_start_pos_in_genome + 2*ERROR, last_successful_location);
					continue;
				}					
				score = c;   // Reseting of score
				edit_dist = ERROR+1;
				VP = ~0; VN = 0;					
				B[0] = B_F[0]; B[1] = B_F[1]; B[2] = B_F[2];B[3] = B_F[3]; B[4] = B_F[4]; 	
				for(k = 0; k < band_len; k++) // verifying for n+2e length    	READ_LENGTH + ERROR + ERROR
				{
					B[0] = B[0] >> 1;
					B[1] = B[1] >> 1;
					B[2] = B[2] >> 1;
					B[3] = B[3] >> 1;
					B[4] = B[4] >> 1;
					if(k + c < READ_LENGTH)
					{
						B[RF[k+c]] = B[RF[k+c]] | MASK;
					}
					X = B[genome[k + verif_start_pos_in_genome]] | VN;
					D0 = ((VP + (X & VP)) ^ VP) | X;
					HN = VP & D0;
					HP = VN | ~(VP | D0);
					X = D0 >> 1;
					VN = X & HP;
					VP = HN | ~(X | HP);

					if(k < (READ_LENGTH-c))
					{
						score = score + 1 - ((D0 >> (W-1)) & 1);
					}
					else
					{
						s = constant1 - k;//s = (W-2) - (k - (READ_LENGTH - c + 1));	
						score = score + ((HP >> s) & 1);
						score = score - ((HN >> s) & 1);
					}
					if(score < edit_dist && (k >= (READ_LENGTH-c)))
					{
						edit_dist = score;
						temp_location = k + verif_start_pos_in_genome;
					}
				}
				if(edit_dist <= ERROR)
				{
					last_successful_location = temp_location - (READ_LENGTH + 2*ERROR);
					endpos_for_mapped_reads[gid + no_of_locations] = temp_location+1;			//Adding 1 as the location is zero-based
					genomic_strand_and_ED_for_mapped_reads[gid + no_of_locations] = 128 + edit_dist;
					no_of_locations = no_of_locations + 1;					
				}				
			}							
		}
		else
		{
			x = x-MIN_QG_LEN;
		}				
		
	}
	//Filtration testing
	eb = EXTRA_BASES; unused_eb = 0;
	if(eb > 0)
	{
		x = READ_LENGTH;	 k_mer_using_eb = 0;
		for (i = 0; i < NUM_OF_QGRAMS; i++)
		{
			x--;
			occurences = 0;
			if(RR[x] != 4)
			{
				sa_start_pos = local_F[RR[x]];
				occurences = local_F[RR[x] + 1] - local_F[RR[x]];
				k = 1; x--;count = 1;			
				while(k)
				{	
					if(RR[x] == 4)			
					{
						if(count < MIN_QG_LEN)
						{
							x = x + count - MIN_QG_LEN;
						}
						occurences = 0;			
						break;					
					}
					temp_location = tally[(sa_start_pos -1)*NO_OF_ALPHABETS + RR[x]];
					occurences = tally[(sa_start_pos + occurences - 1)*NO_OF_ALPHABETS + RR[x]] - temp_location;
					sa_start_pos = local_F[RR[x]] + temp_location;
					if(count < MIN_QG_LEN)
					{
						x--;
						count++;
					}
					else if(occurences <= 1000)
					{
						k=0;
					}
					else if(occurences > 1000)
					{
						if(eb > 0 && RR[x-1] != 4)
						{
							x--;
							eb--;
							k_mer_using_eb = k_mer_using_eb | (1<<i);
						}
						else
						{
							k=0;
						}
					}
				}							
			}
			else
			{
				x = x-MIN_QG_LEN;
			}
		}
		unused_eb = eb; k_mer_using_eb = ~k_mer_using_eb;
	}

	x = READ_LENGTH;	eb = EXTRA_BASES;	// last_successful_location = 0;//no_of_locations=0; 
	for (i = 0; i < NUM_OF_QGRAMS; i++)
	{
		x--;
		occurences = 0;
		if(RR[x] != 4)
		{
			sa_start_pos = local_F[RR[x]];
			occurences = local_F[RR[x] + 1] - local_F[RR[x]];	
			k = 1; x--;count = 1;
			while(k)
			{	
				if(RR[x] == 4)			
				{
					if(count < MIN_QG_LEN)
					{
						x = x + count - MIN_QG_LEN;
					}
					occurences = 0;
					break;					
				}
				temp_location = tally[(sa_start_pos -1)*NO_OF_ALPHABETS + RR[x]];
				occurences = tally[(sa_start_pos + occurences - 1)*NO_OF_ALPHABETS + RR[x]] - temp_location;
				sa_start_pos = local_F[RR[x]] + temp_location;
				if(count < MIN_QG_LEN)
				{
					x--;
					count++;
				}
				else if(occurences <= 1000)
				{
					if((k_mer_using_eb & (1 << i)) > 0 && (unused_eb > 0) && (eb > 0) && (RR[x-1] != 4))
					{
						unused_eb--;
						x--;
						eb--;
						k_mer_using_eb = k_mer_using_eb & (~(1 << i));
					}
					else
					{
						k=0;
					}	
				}
				else if(occurences > 1000)
				{
					if(eb > 0 && RR[x-1] != 4)
					{
						x--;
						eb--;
					}
					else
					{
						k=0;
					}
				}
			}
			occurences = (occurences > 1000)?1000:occurences;
			// no_of_locations = no_of_locations + occurences;
			for(j = 0; j < occurences; j++)
			{
				verif_start_pos_in_genome = SA[sa_start_pos + j] - x - ERROR; //overflow problem can occur if genome length is greater than 2^31 as this variable is signed integer. Keep in mind when dealing with whole genome		
				if(no_of_locations >= CAND_LOC_PER_READ || (verif_start_pos_in_genome > last_successful_location && verif_start_pos_in_genome < last_successful_location  + 4*ERROR))
				{
					// printf("%u 	 %u	%u 		%u	%u 		%u\n",i, x, SA[sa_start_pos + j], verif_start_pos_in_genome, verif_start_pos_in_genome + 2*ERROR, last_successful_location);	
					continue;
				}					
				score = c;   // Reseting of score
				edit_dist = ERROR+1;
				VP = ~0; VN = 0;
				B[0] = B_R[0]; B[1] = B_R[1]; B[2] = B_R[2];B[3] = B_R[3]; B[4] = B_R[4];
				for(k = 0; k < band_len; k++) // verifying for n+2e length    	READ_LENGTH + ERROR + ERROR
				{
					B[0] = B[0] >> 1;
					B[1] = B[1] >> 1;
					B[2] = B[2] >> 1;
					B[3] = B[3] >> 1;
					B[4] = B[4] >> 1;
					if(k + c < READ_LENGTH)
					{
						B[RR[k+c]] = B[RR[k+c]] | MASK;
					}
					X = B[genome[k + verif_start_pos_in_genome]] | VN;
					D0 = ((VP + (X & VP)) ^ VP) | X;
					HN = VP & D0;
					HP = VN | ~(VP | D0);
					X = D0 >> 1;
					VN = X & HP;
					VP = HN | ~(X | HP);

					if(k < (READ_LENGTH-c))
					{
						score = score + 1 - ((D0 >> (W-1)) & 1);
					}
					else
					{
						s = constant1 - k;//s = (W-2) - (k - (READ_LENGTH - c + 1));	
						score = score + ((HP >> s) & 1);
						score = score - ((HN >> s) & 1);
					}
					if(score < edit_dist && (k >= (READ_LENGTH-c)))
					{
						edit_dist = score;
						temp_location = k + verif_start_pos_in_genome;
					}
				}					
				if(edit_dist <= ERROR)
				{
					last_successful_location = temp_location - (READ_LENGTH + 2*ERROR);
					endpos_for_mapped_reads[gid + no_of_locations] = temp_location+1;
					genomic_strand_and_ED_for_mapped_reads[gid + no_of_locations] = edit_dist;
					no_of_locations = no_of_locations + 1;
				}				
			}					
		}
		else
		{
			x = x-MIN_QG_LEN;
		}	
	}//FILTRATION ENDS	
	cand_locs_per_read[get_global_id(0)] = no_of_locations;
	barrier(CLK_GLOBAL_MEM_FENCE);
}


	// if(no_of_locations >= 12000 )//  
	// {
	// 	no_of_locations = 1;
	// }
	// else
	// {
	// 	no_of_locations = 0;	
	// }