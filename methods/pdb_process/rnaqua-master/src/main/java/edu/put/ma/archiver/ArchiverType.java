package edu.put.ma.archiver;

import lombok.Getter;

public enum ArchiverType {

    TAR_GZ(".tar.gz", "tar.gz"), ZIP(".zip", "zip");

    @Getter
    private final String postfix;

    @Getter
    private final String extension;

    ArchiverType(final String postfix, final String extension) {
        this.postfix = postfix;
        this.extension = extension;
    }
}
