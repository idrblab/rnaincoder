package edu.put.ma.archiver;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.rauschig.jarchivelib.ArchiveFormat;
import org.rauschig.jarchivelib.Archiver;
import org.rauschig.jarchivelib.ArchiverFactory;
import org.rauschig.jarchivelib.CompressionType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.io.Files;

public class ArchiverImpl implements edu.put.ma.archiver.Archiver {

    private static final int MINIMAL_ACCEPTABLE_DESTINATION_NAME_SIZE = 3;

    private static final Logger LOGGER = LoggerFactory.getLogger(ArchiverImpl.class);

    private final Archiver archiver;

    ArchiverImpl(final ArchiveFormat format, final CompressionType type) {
        this.archiver = ArchiverFactory.createArchiver(format, type);
    }

    ArchiverImpl(final ArchiveFormat format) {
        this.archiver = ArchiverFactory.createArchiver(format);
    }

    @Override
    public void create(final String archivePath, final File source) {
        try {
            final String archiveName = FilenameUtils.getName(archivePath);
            final File destination = FileUtils.getFile(FilenameUtils.getFullPath(archivePath));
            if (StringUtils.length(destination.getName()) < MINIMAL_ACCEPTABLE_DESTINATION_NAME_SIZE) {
                createWhenDestinationNameIsTooShort(archivePath, source, archiveName);
            } else {
                archiver.create(archiveName, destination, source);
            }
        } catch (IOException e) {
            LOGGER.error(e.getMessage(), e);
        }
    }

    private void createWhenDestinationNameIsTooShort(final String archivePath, final File source,
            final String archiveName) throws IOException {
        final File archiveTempDir = Files.createTempDir();
        try {
            archiver.create(archiveName, archiveTempDir, source);
            Files.move(FileUtils.getFile(archiveTempDir, archiveName), new File(archivePath));
        } finally {
            edu.put.ma.utils.FileUtils.deleteDirectory(archiveTempDir);
        }
    }

    @Override
    public void extract(final File archive, final File destination) {
        try {
            archiver.extract(archive, destination);
        } catch (IOException e) {
            LOGGER.error(e.getMessage(), e);
        }
    }

}
