//CS 576 H4
//Author: Tyler Officer

import java.util.*;
import java.io.*;

public class exon_accuracy {
	
	public static void main(String[] args) throws IOException {
		String truth_file = args[0];
		String prediction_file = args[1];

		ArrayList<String> truth = getSequences(truth_file);
		ArrayList<String> prediction = getSequences(prediction_file);

		getStats(truth, prediction);
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

	public static void getStats(ArrayList<String> truth, ArrayList<String> prediction) {
		//stats[0] = total num positions, stats[1] = total correct predictions
		//stats[2] = correct exonic prediction, stats[3] = true exonic positions, stats[4] = predicted exonic positions
		double[] stats = new double[5];
		char[] t = null;
		char[] p = null;
		int n = truth.size(); //assuming truth length and prediction length equal

		for (int i = 0; i < n; i++) {
			t = truth.get(i).toCharArray();
			p = prediction.get(i).toCharArray();
			
			//similary assuming length of words to be compared is equal
			int m = t.length;
			stats[0] += m;
			for (int j = 0; j < m; j++) {
				if (t[j] == p[j]) stats[1]++;
				if (Character.isUpperCase(t[j]) && t[j] == p[j]) stats[2]++;
				if (Character.isUpperCase(t[j])) stats[3]++;
				if (Character.isUpperCase(p[j])) stats[4]++;
			}
		}
		
		double accuracy = stats[1]/stats[0];
		double recall = stats[2]/stats[3];
		double precision = stats[2]/stats[4];

		System.out.println(String.format("%.3g", accuracy));
		System.out.println(String.format("%.3g", recall));
		System.out.println(String.format("%.3g", precision));
	}
}