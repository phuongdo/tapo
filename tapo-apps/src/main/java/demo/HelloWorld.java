package demo;

import org.cnrs.crbm.lib.io.DataIO;

import java.util.List;

/**
 * Created by pdoviet on 1/27/2015.
 */
public class HelloWorld {

    public static void main(String[] args) {


        // permute(java.util.Arrays.asList(1, 2, 3, 4, 5, 6, 7,8,9), 0);


        List<String> rows = DataIO.readLines("data/CLASSPATH");

        for (String row : rows) {
            System.out.print(row + ";");
        }


//        try {
//
//            String cmd = "/home/pdoviet/save/BioApps/MUSTANG_v3.2.2/bin/mustang-3.2.2 -f /home/pdoviet/work/tmp/input/pdbs.in  -o /home/pdoviet/work/result.out -F FASTA -s OFF";
//            Runtime runtime = Runtime.getRuntime();
//            Process process = runtime.exec(new String[] { "/bin/bash",
//                    "-c", cmd });
//
//            int exitValue = process.waitFor();
//            // System.out.println(cmd);
//            BufferedReader buf = new BufferedReader(new InputStreamReader(
//                    process.getInputStream()));
//            String line = "";
//            while ((line = buf.readLine()) != null) {
//                // builder.append(line + "<br/>");
//                System.out.println("exec response: " + line);
//            }
//        } catch (Exception e) {
//            System.out.println(e);
//        }
    }

    static void permute(java.util.List<Integer> arr, int k) {
        for (int i = k; i < arr.size(); i++) {
            java.util.Collections.swap(arr, i, k);
            permute(arr, k + 1);
            java.util.Collections.swap(arr, k, i);
        }
        if (k == arr.size() - 1) {

            double total = arr.get(0) + ((double) 13 * arr.get(1) / arr.get(2)) + arr.get(3) + ((double) 12 * arr.get(4)) - arr.get(5) + (arr.get(6) * arr.get(7) / arr.get(8));

            if (total == 88.00) {
                System.out.println(java.util.Arrays.toString(arr.toArray()));
            }

        }
    }

}
