package edu.put.ma.utils;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.io.FileUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.io.Files;

import edu.put.ma.model.services.PdbStructure;
import edu.put.ma.model.services.StructuresSet;

public final class PdbReader {

    private static final String PDB_EXTENSION = ".pdb";

    private static final Logger LOGGER = LoggerFactory.getLogger(PdbReader.class);

    private PdbReader() {
        // hidden constructor
    }

    public static final void readSingleStructure(final String singleStructurePath,
            final StructuresSet structuresSet, final boolean isTarget) {
        final Map<Integer, List<String>> atoms = new HashMap<Integer, List<String>>();
        final File structureFile = FileUtils.getFile(singleStructurePath);
        final int modelNo = readStructure(atoms, structureFile);
        addStructure(structuresSet, atoms, modelNo, structureFile, isTarget);
    }

    public static final void readMultiplePdbModelsDirectory(final String multiplePdbModelsDirectoryPath,
            final StructuresSet structuresSet) {
        final File dir = FileUtils.getFile(multiplePdbModelsDirectoryPath);
        final File[] files = dir.listFiles();
        for (int i = 0; i < files.length; i++) {
            if (files[i].getName().endsWith(PDB_EXTENSION)) {
                readSingleStructure(files[i].getAbsolutePath(), structuresSet, false);
            }
        }
    }

    private static int readStructure(final Map<Integer, List<String>> atoms, final File structureFile) {
        int modelNo = 1;
        try {
            for (String line : Files.readLines(structureFile, Charset.defaultCharset())) {
                if (line.startsWith("MODEL")) {
                    modelNo = Integer.parseInt(org.apache.commons.lang3.StringUtils
                            .trim(org.apache.commons.lang3.StringUtils.substring(line, 10, 14)));
                } else if ((line.startsWith("ATOM")) || (line.startsWith("HETATM"))) {
                    initContainerForParticularModel(atoms, modelNo);
                    atoms.get(modelNo).add(line);
                }
            }
        } catch (IOException e) {
            LOGGER.error(e.getMessage(), e);
        }
        return modelNo;
    }

    private static void initContainerForParticularModel(final Map<Integer, List<String>> atoms,
            final int modelNo) {
        if (!atoms.containsKey(modelNo)) {
            atoms.put(modelNo, new ArrayList<String>());
        }
    }

    private static void addStructure(final StructuresSet structuresSet,
            final Map<Integer, List<String>> atoms, final int modelNo, final File structureFile,
            final boolean isTarget) {
        if (CollectionUtils.size(atoms) == 1) {
            if (isTarget) {
                structuresSet.addStructure(getPdbStructure(atoms, modelNo, structureFile.getName()), 0);
            } else {
                structuresSet.addStructure(getPdbStructure(atoms, modelNo, structureFile.getName()));
            }
        } else {
            final SortedSet<Integer> keys = new TreeSet<Integer>(atoms.keySet());
            int index = 0;
            for (Integer key : keys) {
                if (isTarget) {
                    structuresSet.addStructure(
                            getPdbStructure(
                                    atoms,
                                    key,
                                    structureFile.getName().replaceAll("\\.pdb",
                                            "#" + key.toString() + PDB_EXTENSION)), index++);
                } else {
                    structuresSet.addStructure(getPdbStructure(atoms, key, structureFile.getName()
                            .replaceAll("\\.pdb", "#" + key.toString() + PDB_EXTENSION)));
                }
            }
        }
    }

    private static PdbStructure getPdbStructure(final Map<Integer, List<String>> atoms, final int modelNo,
            final String filename) {
        final PdbStructure structure = new PdbStructure();
        structure.setFilename(filename);
        structure.setAtomList(atoms.get(modelNo));
        return structure;
    }

}
