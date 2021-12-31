package org.cnrs.crbm.lib.utils;

import org.cnrs.crbm.lib.conf.ConfigUtil;
import org.cnrs.crbm.lib.conf.Path;

import java.io.BufferedReader;
import java.io.InputStreamReader;

public class ExecUtils {

    private static final String BASH_EXECUTABLE = Path.CYGWIN_PATH
            + ConfigUtil.getInstance().getProperty("BASH_EXECUTABLE");

    public static String execShellCmdLinux(String cmd) {
        // StringBuilder builder = new StringBuilder();
        StringBuffer buffer = new StringBuffer();
        try {
            Runtime runtime = Runtime.getRuntime();
            Process process = runtime.exec(new String[]{BASH_EXECUTABLE,
                    "-c", cmd});

            int exitValue = process.waitFor();
            // System.out.println(cmd);

            BufferedReader buf = new BufferedReader(new InputStreamReader(
                    process.getInputStream()));
            String line = "";
            while ((line = buf.readLine()) != null) {
                buffer.append(line + "\n");
                // System.out.println("exec response: " + line);
            }
            process.destroy();
        } catch (Exception e) {
            //System.out.println(e);
        }
        return buffer.toString();
    }

}
