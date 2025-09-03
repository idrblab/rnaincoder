package edu.put.ma.archiver;

import java.io.File;

import org.apache.commons.lang3.StringUtils;
import org.rauschig.jarchivelib.ArchiveFormat;
import org.rauschig.jarchivelib.CompressionType;

public final class ArchiverFactory {

    private ArchiverFactory() {
        // hidden constructor
    }

    public static final Archiver getArchiver(final ArchiverType archiverType) {
        return (archiverType == ArchiverType.TAR_GZ) ? new ArchiverImpl(ArchiveFormat.TAR,
                CompressionType.GZIP) : new ArchiverImpl(ArchiveFormat.ZIP);
    }

    public static final Archiver getArchiver(final File file) {
        return getArchiver(file.getName());
    }

    public static final Archiver getArchiver(final String fileName) {
        return (StringUtils.endsWithIgnoreCase(fileName, ArchiverType.TAR_GZ.getPostfix())) ? new ArchiverImpl(
                ArchiveFormat.TAR, CompressionType.GZIP) : new ArchiverImpl(ArchiveFormat.ZIP);
    }

    public static final boolean isArchive(final File file) {
        return isArchive(file.getName());
    }

    public static final boolean isArchive(final String fileName) {
        return StringUtils.endsWithIgnoreCase(fileName, ArchiverType.ZIP.getPostfix())
                || StringUtils.endsWithIgnoreCase(fileName, ArchiverType.TAR_GZ.getPostfix());
    }

    public static final String getArchiverPostfix(final File file) {
        return getArchiverPostfix(file.getName());
    }

    public static final String getArchiverPostfix(final String fileName) {
        return StringUtils.endsWithIgnoreCase(fileName, ArchiverType.TAR_GZ.getPostfix()) ? ArchiverType.TAR_GZ
                .getPostfix() : ArchiverType.ZIP.getPostfix();
    }

}
