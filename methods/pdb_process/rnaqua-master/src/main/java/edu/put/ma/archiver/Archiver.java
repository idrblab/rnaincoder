package edu.put.ma.archiver;

import java.io.File;

public interface Archiver {

    void create(String archivePath, File source);

    void extract(File archive, File destination);
}
