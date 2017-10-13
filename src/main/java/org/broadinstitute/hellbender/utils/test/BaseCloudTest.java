package org.broadinstitute.hellbender.utils.test;

import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

/**
 * Base class for tests that use cloud resources, we include
 * a few helpful methods.
 */
public abstract class BaseCloudTest extends BaseTest {

  protected FileSystem getAuthenticatedGcs(String bucket) throws IOException {
    byte[] creds = Files.readAllBytes(Paths.get(getGoogleServiceAccountKeyPath()));
    return BucketUtils.getAuthenticatedGcs(getGCPTestProject(), bucket, creds);
  }

  protected void helpDebugAuthError() {
    final String key = "GOOGLE_APPLICATION_CREDENTIALS";
    String credsFile = System.getenv(key);
    if (null == credsFile) {
      System.err.println("$"+key+" is not defined.");
      return;
    }
    System.err.println("$"+key+" = " + credsFile);
    Path credsPath = Paths.get(credsFile);
    boolean exists = Files.exists(credsPath);
    System.err.println("File exists: " + exists);
    if (exists) {
      try {
        System.err.println("Key lines from file:");
        printKeyLines(credsPath, "\"type\"", "\"project_id\"", "\"client_email\"");
      } catch (IOException x2) {
        System.err.println("Unable to read: " + x2.getMessage());
      }
    }
  }

  private void printKeyLines(Path path, String... keywords) throws IOException {
    for (String line : Files.readAllLines(path)) {
      for (String keyword : keywords) {
        if (line.contains(keyword)) {
          System.err.println(line);
        }
      }
    }
  }

}
