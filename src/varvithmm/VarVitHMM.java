/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package varvithmm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

/**
 *
 * @author lucas
 */
public class VarVitHMM {
    final static String DESCRIPTION = "Program requires 3 parameters and can be run as follows:\n" +
"   java -jar varVitHMM.jar [FILE1] [FILE2] [FILE3]\n\n" +
" EX.\n" +
"    java -jar varVitHMM.jar hmm-tm.txt sequences-project2.txt sequences-project2-viterbi-output.txt\n\n" +
"   [FILE1] - file with HMM parameters\n" +
"   [FILE2] - input file with sequences in fasta format\n" +
"   [FILE3] - name of an output file.\n\n" +

"";
    String[] obs;
    char[] hid;
    HashMap<String, Integer> obsMap = new HashMap<>();
    
    double[] init;
    double[][] trans;
    double[][] em;
    int[] varD;
    
    ArrayList<String> seqs = new ArrayList<>();
    ArrayList<String> names = new ArrayList<>();
    String[] states;
    
    int K, N;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        if(args.length!=3){
            System.out.println(DESCRIPTION);
            System.exit(0);
        }
        else{
            VarVitHMM m = new VarVitHMM();
            m.VarVitHMM(args);
        }
        // TODO code application logic here
    }
    public double log(double x){
        if(x==0){
            return -Double.MAX_VALUE;
        }
        else{
            return Math.log(x);
        }
    }

    public void VarVitHMM(String[] args){
        readMM(args[0]);
        readInp(args[1]);
        K=init.length;
        long startTime = System.currentTimeMillis();
        for(int m=0; m<seqs.size();m++){
            char[] X = seqs.get(m).toCharArray();
            N = X.length;
            long startTimeOmega = System.currentTimeMillis();
            double[][] omega = new double[K][N];
            //fill omega table
            for(int k=0;k<K;k++){
                omega[k][0] = init[k] + em[k][obsMap.get(String.valueOf(X[0]))];
//                omega[k][0] = log(init[k]) + log(em[k][obsMap.get(String.valueOf(X[0]))]);
            }
            for(int k=0;k<K;k++){
                
                
                if(varD[k]==1){
//                    omega[k][1] = getMaxOmega(omega, 0, k) + log(em[k][obsMap.get(String.valueOf(X[1]))]);
                    omega[k][1] = getMaxOmega(omega, 0, k) + em[k][obsMap.get(String.valueOf(X[1]))];
                }
                else{
                    omega[k][1]=omega[k][0];
                }
            }
            for(int k=0;k<K;k++){
                if(varD[k]==1){
                    omega[k][2] = getMaxOmega(omega, 1, k) + em[k][obsMap.get(String.valueOf(X[2]))];
//                    omega[k][2] = getMaxOmega(omega, 1, k) + log(em[k][obsMap.get(String.valueOf(X[2]))]);
                }
                else if(varD[k]==2){
//                    omega[k][2] = getMaxOmega(omega, 0, k) + log(em[k][obsMap.get(String.valueOf(X[1])+String.valueOf(X[2]))]);                    
                    omega[k][2] = getMaxOmega(omega, 0, k) + em[k][obsMap.get(String.valueOf(X[1])+String.valueOf(X[2]))];                    
                }
                else{
                    omega[k][2]=omega[k][1];
                }
            }
            
            
            
            for(int n=3;n<N;n++){
                for(int k=0;k<K;k++){
                    if(varD[k]==1){
                        omega[k][n] = getMaxOmega(omega, n-1, k) + em[k][obsMap.get(String.valueOf(X[n]))];
//                        omega[k][n] = getMaxOmega(omega, n-1, k) + log(em[k][obsMap.get(String.valueOf(X[n]))]);
                    }
                    else{
                        String symbol="";
                        for(int d=varD[k]-1; d>=0;d--){
                            symbol+=String.valueOf(X[n-d]);
                        }
                        omega[k][n]=getMaxOmega(omega, n-varD[k], k)+em[k][obsMap.get(symbol)];
//                        omega[k][n]=getMaxOmega(omega, n-varD[k], k)+log(em[k][obsMap.get(symbol)]);
                    }
//                    if(varD[k]==1){
//                        omega[k][n] = getMaxOmega(omega, n-1, k) + log(em[k][obsMap.get(String.valueOf(X[n]))]);
//                    }
//                    else if(varD[k]==2){
//                        omega[k][n] = getMaxOmega(omega, n-2, k) + log(em[k][obsMap.get(String.valueOf(X[n-1])+String.valueOf(X[n]))]);                    
//                    }
//                    else{
//                        omega[k][n] = getMaxOmega(omega, n-3, k) + log(em[k][obsMap.get(String.valueOf(X[n-2])+String.valueOf(X[n-1])+String.valueOf(X[n]))]);                    
//                    }
                }
            }
            long timeOmega = System.currentTimeMillis()-startTimeOmega;
            
            //backtrack
            long startTimeBack = System.nanoTime();
            char[] Z = new char[N];
            //find last hidden state
            int nextState = getMaxOmegaIdx(omega);
            double omegaLogJP = omega[nextState][N-1];
            int end = N-2;
            if(varD[nextState]==1){ //if last state emitts 1 symbol
                Z[N-1] = hid[nextState];
            }
            if(varD[nextState]==2){ //if last state emitt 2 symbols
                Z[N-1] = hid[nextState];
                Z[N-2] = hid[nextState];
                end=N-3;
            }
            else{ //if last state emitt 3 symbols
                Z[N-1] = hid[nextState];
                Z[N-2] = hid[nextState];
                Z[N-3] = hid[nextState];
                end=N-4;                
            }
            //find the rest of hidden states
            for(int n=end;n>=2;n--){
                nextState = getMaxIdx(omega, nextState, n);
                if(varD[nextState]==1){
                    Z[n] = hid[nextState];
                }
                else if(varD[nextState]==2){
                    Z[n] = hid[nextState];
                    Z[n-1] = hid[nextState];
                    n-=1;
                }
                else{
                    Z[n] = hid[nextState];
                    Z[n-1] = hid[nextState];
                    Z[n-2] = hid[nextState];
                    n-=2;                
                }
            }
            //find first 2 hidden states
            nextState = getMaxIdx(omega, nextState, 1);
            if(varD[nextState]==1){
                Z[1] = hid[nextState];
                nextState = getMaxIdx(omega, nextState, 0);
                if(varD[nextState]==1){
                    Z[0] = hid[nextState];
                }
            }
            else if(varD[nextState]==2){
                Z[1] = hid[nextState];
                Z[0] = hid[nextState];
            }
            float timeBack = Float.valueOf(System.nanoTime()-startTimeBack)/1000000;
            //saveTimings(timeOmega, timeBack, N, K);
            states[m] = String.valueOf(Z);
//            double logJP = getLogJP(X,Z);
//            if(omegaLogJP == logJP){
//                System.out.println("OK");
//            } 
//            else{
//                System.out.println(omegaLogJP + " "+ logJP);
//            }
        }
        long runTime = System.currentTimeMillis()-startTime;
        System.out.println(states.length+" sequences decoded in "+(double)runTime/1000+"s.");
        saveOut(args[2]);
    }
    
    //saves running time for timing experiment
    private void saveTimings(long omegaTime, float backTime, int N, int K){
        try {
            File file = new File("timingsKrandom.csv");
            FileWriter fw = new FileWriter(file.getAbsoluteFile(),true);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write(N+", "+K+", "+omegaTime+", "+backTime);
            bw.newLine();
            bw.close();
        } catch (IOException e) {
                e.printStackTrace();
        }               
        
    }
    
//    private double getLogJP(char[] X, char[] Z){
//        double logJP = log(init[hidMap.get(Z[0])]) + log(em[hidMap.get(Z[0])][obsMap.get(String.valueOf(X[0]))]);
//        for(int n=1;n<N;n++){
//            logJP+= log(trans[hidMap.get(Z[n-1])][hidMap.get(Z[n])]);
//            logJP+= log(em[hidMap.get(Z[n])][obsMap.get(String.valueOf(X[n]))]);
//        }
//        return logJP;
//    }
    
    //returns max value of previous column + transition to current cell
    private double getMaxOmega(double[][] omega, int n, int k){
        double max = Double.NEGATIVE_INFINITY;
        for(int j=0;j<K;j++){
            double d = omega[j][n] + trans[j][k];
//            double d = omega[j][n] + log(trans[j][k]);
            if(d>max){
                max=d;
            }
        }
        return max;
    }
    
    private int getMaxIdx(double[][] omega, int nextState, int n){
        int idx = -1;
        double max = Double.NEGATIVE_INFINITY;
        for(int k=0;k<K;k++){
            double d = omega[k][n] + trans[k][nextState];
//            double d = omega[k][n] + log(trans[k][nextState]);
            if(max<d){
                max = d;
                idx = k;
            }
        }
        return idx;
    }
    
    //returns max row index of last column,
    //most probable state at the end of seq.
    private int getMaxOmegaIdx(double[][] omega){
        double max = Double.NEGATIVE_INFINITY;
        int idx = -1;
        for(int k=0;k<K;k++){
            if(max < omega[k][omega[k].length-1]){
                max = omega[k][omega[k].length-1];
                idx = k;
            }
        }
        return idx;
    }    
    
    private void readMM(String input){
        try{
            BufferedReader br1 = new BufferedReader(new FileReader(input));
            String line = br1.readLine();
            while(line!=null){
                if(line.contains("Init")){
                    String[] l = br1.readLine().split(" ");
                    init = new double[l.length];
                    for(int i=0;i<init.length;i++){
                        init[i] = log(Double.valueOf(l[i]));
                    }
                }
                if(line.contains("Transitions")){
                    String[] l = br1.readLine().split(" ");
                    int n = l.length;
                    trans = new double[n][n];
                    for(int i=0;i<n;i++){
                        for(int j=0;j<n;j++){
                            trans[i][j] = log(Double.valueOf(l[j]));
                        }
                        l = br1.readLine().split(" ");
                    }
                }
                if(line.contains("Emissions")){
                    String[] l = br1.readLine().split(" ");
                    int n = l.length;
                    em = new double[init.length][n];
                    for(int i=0;i<init.length;i++){
                        for(int j=0;j<n;j++){
                            em[i][j] = log(Double.valueOf(l[j]));
                        }
                        line = br1.readLine();
                        if(line!=null)
                            l = line.split(" ");
                        
                    }
                }
                else if(line.contains("hidden")){
                    String[] l = br1.readLine().split(" ");
                    int n = l.length;
                    hid = new char[n];
                    for(int i=0;i<n;i++){
                        hid[i]=l[i].charAt(0);
                    }                    
                }
                else if(line.contains("observables")){
                    String[] l = br1.readLine().split(" ");
                    int n = l.length;
                    obs = new String[n];
                    for(int i=0;i<n;i++){
                        obs[i]=l[i];
                        obsMap.put(l[i], i);
                    }                    
                }
                else if(line.contains("variable")){
                    String[] l = br1.readLine().split(" ");
                    int n = l.length;
                    varD = new int[n];
                    for(int i=0;i<n;i++){
                        varD[i]=Integer.valueOf(l[i]);
                    }                    
                }
                line = br1.readLine();
            }
        }catch(Exception e){
            e.printStackTrace();
        }                
    }

    private void readInp(String input){
        System.out.println("Reading input sequences from: "+ input);

        try{
            BufferedReader br1 = new BufferedReader(new FileReader(input));
            String line = br1.readLine();
            seqs = new ArrayList<>();
            while(line!=null){
                if(line.startsWith(">")){
                    String[] l = line.split(" ");
                    names.add(l[0]);
                    line = br1.readLine();
                }
                else{
                    String l = line;
                    line = br1.readLine();
                    while(line!=null && !line.startsWith(">")){
                        l+=line;
                        line = br1.readLine();
                    }
                    while(l.indexOf('N')!=-1){
                        l=l.replaceFirst("N", randomChar('N'));
                    }
                    while(l.indexOf('Y')!=-1){
                        l=l.replaceFirst("Y", randomChar('Y'));
                    }
                    while(l.indexOf('K')!=-1){
                        l=l.replaceFirst("K", randomChar('K'));
                    }
                    while(l.indexOf('R')!=-1){
                        l=l.replaceFirst("R", randomChar('R'));
                    }
                    while(l.indexOf('M')!=-1){
                        l=l.replaceFirst("M", randomChar('M'));
                    }
                    while(l.indexOf('S')!=-1){
                        l=l.replaceFirst("S", randomChar('S'));
                    }
                    while(l.indexOf('W')!=-1){
                        l=l.replaceFirst("W", randomChar('W'));
                    }
                    while(l.indexOf('B')!=-1){
                        l=l.replaceFirst("B", randomChar('B'));
                    }
                    while(l.indexOf('D')!=-1){
                        l=l.replaceFirst("D", randomChar('D'));
                    }
                    while(l.indexOf('H')!=-1){
                        l=l.replaceFirst("H", randomChar('H'));
                    }
                    while(l.indexOf('V')!=-1){
                        l=l.replaceFirst("V", randomChar('V'));
                    }
                    seqs.add(l);
                }
            }
        }catch(Exception e){
            e.printStackTrace();
        }
        states = new String[seqs.size()];
    }
    
    private String randomChar(char c){
        Random r = new Random();
        String n="";
        if(c=='N'){
            n = String.valueOf(obs[r.nextInt(4)]);
        }
        if(c=='Y'){
            char[] nuc = new char[]{'C', 'T'};
            n = String.valueOf(nuc[r.nextInt(2)]);
        }
        if(c=='K'){
            char[] nuc = new char[]{'G', 'T'};
            n = String.valueOf(nuc[r.nextInt(2)]);
        }
        if(c=='R'){
            char[] nuc = new char[]{'G', 'A'};
            n = String.valueOf(nuc[r.nextInt(2)]);
        }
        if(c=='M'){
            char[] nuc = new char[]{'A', 'C'};
            n = String.valueOf(nuc[r.nextInt(2)]);
        }
        if(c=='S'){
            char[] nuc = new char[]{'G', 'C'};
            n = String.valueOf(nuc[r.nextInt(2)]);
        }
        if(c=='W'){
            char[] nuc = new char[]{'A', 'T'};
            n = String.valueOf(nuc[r.nextInt(2)]);
        }
        if(c=='B'){
            char[] nuc = new char[]{'G', 'T', 'C'};
            n = String.valueOf(nuc[r.nextInt(3)]);
        }
        if(c=='D'){
            char[] nuc = new char[]{'G', 'A', 'T'};
            n = String.valueOf(nuc[r.nextInt(3)]);
        }
        if(c=='H'){
            char[] nuc = new char[]{'A', 'C', 'T'};
            n = String.valueOf(nuc[r.nextInt(3)]);
        }
        if(c=='V'){
            char[] nuc = new char[]{'G', 'C', 'A'};
            n = String.valueOf(nuc[r.nextInt(3)]);
        }
        return n;
    }

    private void saveOut(String output){
        try {
            File file = new File(output);
            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);
            for(int i=0;i<states.length;i++){
                bw.write(names.get(i));
                bw.newLine();
                bw.write(states[i]);
                bw.newLine();
            }
            bw.close();
        } catch (IOException e) {
                e.printStackTrace();
        }               
    }        
   
}
