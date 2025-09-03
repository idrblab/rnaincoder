package edu.put.ma;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.util.Properties;

import javax.ws.rs.client.Entity;
import javax.ws.rs.client.WebTarget;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.base.Preconditions;

import edu.put.ma.model.input.CommonInputModel;
import edu.put.ma.model.input.ComparisonInputModelImpl;
import edu.put.ma.model.input.InputModelsFactory;
import edu.put.ma.model.services.StructuresSet;
import edu.put.ma.templates.Templates;
import edu.put.ma.templates.TemplatesImpl;
import edu.put.ma.utils.PreconditionUtils;
import edu.put.ma.utils.ServicesProvider;

public class Assessment {

    public static final String SUPPORTED_IMAGE_EXTENSION = "svg";

    private static final Logger LOGGER = LoggerFactory.getLogger(Assessment.class);

    private Properties properties;

    private Templates templates = new TemplatesImpl();

    private WebTarget webTarget = null;

    Assessment() {
        this.properties = new Properties();
        loadProperties();
        this.webTarget = ServicesProvider.getWebTarget(getServicesProvider());
        initTemplates();
    }

    private void initTemplates() {
        if (!templates.isInitiated()) {
            templates.load(this.getTemplates());
        }
    }

    private String getTemplates() {
        return webTarget.path("services").path("rnalyzer").path("templates").request(MediaType.TEXT_PLAIN)
                .get(String.class);
    }

    private void loadProperties() {
        InputStream inputStream = null;
        try {
            inputStream = this.getClass().getResourceAsStream("app.properties");
            properties.load(inputStream);
        } catch (IOException e) {
            LOGGER.error(e.getMessage(), e);
        } finally {
            IOUtils.closeQuietly(inputStream);
        }
    }

    private String getPropertyValueByKey(final String key) {
        final String propertyValue = properties.getProperty(key);
        PreconditionUtils.checkIfStringIsBlank(propertyValue, String.format("Value of property [%s]", key));
        return propertyValue;
    }

    private String getArtifactId() {
        return getPropertyValueByKey("artifactId");
    }

    private String getServicesProvider() {
        return getPropertyValueByKey("servicesProvider");
    }

    public static void main(final String[] args) {
        final Assessment app = new Assessment();
        app.perform(args);
    }

    public void perform(final String[] args) {
        final CommonInputModel inputModel = InputModelsFactory.getInputModel(args, getArtifactId());
        Preconditions.checkNotNull(inputModel, "Input model");
        if (inputModel.isInputInitializedProperly()) {
            LOGGER.info(inputModel.getInputModelString());
            final String outputFileExtension = StringUtils.lowerCase(FilenameUtils.getExtension(inputModel
                    .getOutputFilePath()));
            final String serviceCode = getServiceCode(inputModel, outputFileExtension,
                    inputModel instanceof ComparisonInputModelImpl);
            process(inputModel.getCommand(), inputModel.getStructuresSet(), serviceCode,
                    inputModel.getOutputFilePath(), outputFileExtension);
        } else {
            inputModel.printHelp(getArtifactId());
        }
    }

    private void process(final Command command, final StructuresSet structuresSet, final String serviceCode,
            final String outputFilePath, final String outputFileExtension) {
        if (structuresSet.getStructureListSize() > 0) {
            Response response = webTarget.path("services").path("rnalyzer").request().post(Entity.json(null));
            final URI location = response.getLocation();
            final String resourceUri = StringUtils.difference(webTarget.getUri().toString(),
                    location.toString());
            try {
                String responseString = "Service does not exist or server is shutdown!";
                if (location != null) {
                    response = webTarget.path(resourceUri).path(serviceCode)
                            .request(getMediaType(command, outputFileExtension))
                            .put(Entity.xml(structuresSet));
                }
                if (response.getStatus() == Response.Status.NO_CONTENT.getStatusCode()) {
                    responseString = "Service is turned off. Try to come later!";
                    LOGGER.info(responseString);
                } else if (response.getStatus() == Response.Status.NOT_FOUND.getStatusCode()) {
                    responseString = "Resource you are trying to connect does not exist. Try to initialize it first!";
                    LOGGER.info(responseString);
                } else {
                    saveOutput(command, outputFilePath, response);
                    LOGGER.info("Command processed properly");
                }
            } catch (Exception e) {
                LOGGER.error(e.getMessage(), e);
            } finally {
                webTarget.path(resourceUri).request().delete();
            }
        }
    }

    private void saveOutput(final Command command, final String outputFilePath, final Response response) {
        if ((CommandEnum.Structure3d.ORIGINAL_3D == command)
                || (CommandEnum.Structure3d.RENUMERATED_3D == command)
                || (CommandEnum.DEFORMATION_PROFILE == command)) {
            writeBinaryOutputToFile(outputFilePath, response);
        } else {
            writeTextualOutputToFile(outputFilePath, response.readEntity(String.class));
        }
    }

    private static void writeBinaryOutputToFile(final String outputPath, final Response response) {
        FileOutputStream fos = null;
        try {
            byte[] bytesArray = IOUtils.toByteArray(response.readEntity(InputStream.class));
            fos = new FileOutputStream(FileUtils.getFile(outputPath));
            fos.write(bytesArray);
        } catch (IOException e) {
            LOGGER.error(e.getMessage(), e);
        } finally {
            IOUtils.closeQuietly(fos);
        }
    }

    private static void writeTextualOutputToFile(final String outputFilePath, final String res) {
        try {
            FileUtils.write(FileUtils.getFile(outputFilePath), res);
        } catch (IOException e) {
            LOGGER.error(e.getMessage(), e);
        }
    }

    private static String getMediaType(final Command command, final String extension) {
        if ((CommandEnum.Structure3d.ORIGINAL_3D == command)
                || (CommandEnum.Structure3d.RENUMERATED_3D == command)) {
            return "application/zip";
        } else if (CommandEnum.DEFORMATION_PROFILE == command) {
            if (SUPPORTED_IMAGE_EXTENSION.equals(extension)) {
                return "image/png";
            } else {
                return "application/zip";
            }
        }
        return MediaType.APPLICATION_XML;
    }

    private static String getServiceCode(final CommonInputModel inputModel, final String extension,
            final boolean isComparisonWithTarget) {
        final StringBuilder serviceCodeBuilder = new StringBuilder(inputModel.getServiceCode());
        if ((CommandEnum.DEFORMATION_PROFILE == inputModel.getCommand())
                && (SUPPORTED_IMAGE_EXTENSION.equals(extension))) {
            serviceCodeBuilder.append("i");
        } else if (CommandEnum.getSequenceAndStructure3dServices().contains(inputModel.getCommandName())
                && isComparisonWithTarget) {
            serviceCodeBuilder.append("s");
        }
        return serviceCodeBuilder.toString();
    }
}
