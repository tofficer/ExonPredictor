//Bioinformatics
//Tyler Officer

import java.util.*;
import java.io.*;

public class predict_exons {
	
	public static void main(String[] args) throws IOException {
		String train_file = args[0];
		String test_file = args[1];

		ArrayList<String> train_sequences = getSequences(train_file);
		ArrayList<String> test_sequences = getSequences(test_file);

		double[] hmm_params = train(train_sequences);
		test(test_sequences, hmm_params);
	}

	public static double[] train(ArrayList<String> sequences) {
		//use laplace smoothing ie start counts at 1
		double[] hmm_params = new double[]{1,1,1,1,1,1,1,1,1,1,1,1};

		for (String seq : sequences) {
			boolean prev_uppercase = true;
			boolean curr_uppercase = true;
			for (int i = 0; i < seq.length(); i++) {
				char c = seq.charAt(i);
				switch(c) {
					case 'A': hmm_params[4]++; curr_uppercase = true; break;
					case 'C': hmm_params[5]++; curr_uppercase = true; break;
					case 'G': hmm_params[6]++; curr_uppercase = true; break;
					case 'T': hmm_params[7]++; curr_uppercase = true; break;
					case 'a': hmm_params[8]++; curr_uppercase = false; break;
					case 'c': hmm_params[9]++; curr_uppercase = false; break;
					case 'g': hmm_params[10]++; curr_uppercase = false; break;
					case 't': hmm_params[11]++; curr_uppercase = false; break;
				}

				if (i == 0) continue;

				if (prev_uppercase && curr_uppercase) hmm_params[0]++;
				else if (prev_uppercase) hmm_params[1]++;
				else if (curr_uppercase) hmm_params[2]++;
				else hmm_params[3]++;

				prev_uppercase = curr_uppercase;
			}
		}
		
		double ex_trans = hmm_params[0] + hmm_params[1] + sequences.size() + 1;
		double in_trans = hmm_params[2] + hmm_params[3];
		double ex_emit = hmm_params[4] + hmm_params[5] + hmm_params[6] + hmm_params[7];
		double in_emit = hmm_params[8] + hmm_params[9] + hmm_params[10] + hmm_params[11];

		//from start transitions
		//being to exon has 1 prob @ beginning of character sequence

		//from exon transitions
		hmm_params[0] = Math.log(hmm_params[0]/ex_trans); //exon to exon
		hmm_params[1] = Math.log(hmm_params[1]/ex_trans); //exon to intron
		//exon to end has 1 prob @ end of char seq

		//from intron transitions
		hmm_params[2] = Math.log(hmm_params[2]/in_trans); //intron to exon
		hmm_params[3] = Math.log(hmm_params[3]/in_trans); //intron to intron
		
		//exon emissions
		hmm_params[4] = Math.log(hmm_params[4]/ex_emit); //A
		hmm_params[5] = Math.log(hmm_params[5]/ex_emit); //C
		hmm_params[6] = Math.log(hmm_params[6]/ex_emit); //G
		hmm_params[7] = Math.log(hmm_params[7]/ex_emit); //T
		
		//intron emissions
		hmm_params[8] = Math.log(hmm_params[8]/in_emit); //a
		hmm_params[9] = Math.log(hmm_params[9]/in_emit); //c
		hmm_params[10] = Math.log(hmm_params[10]/in_emit); //g
		hmm_params[11] = Math.log(hmm_params[11]/in_emit); //t

		return hmm_params;
	}

	//using viterbi algo
	public static void test(ArrayList<String> sequences, double[] hmm_params) {
		double[][] viterbi = null;
		int[][] ptrs = null;

		for (String seq : sequences) {
			int n = seq.length();
			
			viterbi = new double[2][n];
			ptrs = new int[2][n]; //keep track of which state you're coming from
			
			for (int i = 0; i < n; i++) {
				char c = seq.charAt(i);
				switch(c) {
					case 'A': viterbi[0][i] = hmm_params[4]; viterbi[1][i] = hmm_params[8]; break;
					case 'C': viterbi[0][i] = hmm_params[5]; viterbi[1][i] = hmm_params[9]; break;
					case 'G': viterbi[0][i] = hmm_params[6]; viterbi[1][i] = hmm_params[10]; break;
					case 'T': viterbi[0][i] = hmm_params[7]; viterbi[1][i] = hmm_params[11]; break;
				}

				if (i < 1) continue;

				//to exon state
				if ( (i==1) || ((viterbi[0][i-1]+hmm_params[0]) > (viterbi[1][i-1]+hmm_params[2])) ) {
					ptrs[0][i] = 0;
					viterbi[0][i] += (viterbi[0][i-1]+hmm_params[0]);
				} else {
					ptrs[0][i] = 1;
					viterbi[0][i] += (viterbi[1][i-1]+hmm_params[2]);
				}

				//to intron state
				if ( (i==1) || ((viterbi[0][i-1]+hmm_params[1]) > (viterbi[1][i-1]+hmm_params[3])) ) {
					ptrs[1][i] = 0;
					viterbi[1][i] += (viterbi[0][i-1]+hmm_params[1]);
				} else {
					ptrs[1][i] = 1;
					viterbi[1][i] += (viterbi[1][i-1]+hmm_params[3]); 
				}

			}

			StringBuilder sb = new StringBuilder();
			int traceback = 0; //there is only a transition from exon to end
			for (int i = n-1; i >= 0; i--) {
				char ch = seq.charAt(i);
				if (traceback == 1) sb.insert(0, Character.toLowerCase(ch));
				else sb.insert(0, ch);

				traceback = ptrs[traceback][i];
			}
			
			System.out.println(sb.toString());
		}
	}

	public static ArrayList<String> getSequences(String infile) throws IOException, FileNotFoundException {
		ArrayList<String> list = new ArrayList<String>();
		try (BufferedReader br = new BufferedReader(new FileReader(infile))) {
			String line;
			while ((line = br.readLine()) != null) list.add(line);
		}
		//br.close();

		return list;
	}

}
