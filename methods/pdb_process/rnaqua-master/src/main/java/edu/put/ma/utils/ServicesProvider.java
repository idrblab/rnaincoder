package edu.put.ma.utils;

import java.security.KeyStore;
import java.security.KeyStoreException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;

import javax.net.ssl.SSLContext;
import javax.net.ssl.TrustManager;
import javax.net.ssl.TrustManagerFactory;
import javax.net.ssl.X509TrustManager;
import javax.ws.rs.client.Client;
import javax.ws.rs.client.ClientBuilder;
import javax.ws.rs.client.WebTarget;
import javax.ws.rs.core.UriBuilder;

import org.glassfish.jersey.client.ClientConfig;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.base.Preconditions;
import com.google.common.collect.Iterables;

public final class ServicesProvider {

    private static final String PUTWSS = "PUTWSs";

    private static final Logger LOGGER = LoggerFactory.getLogger(ServicesProvider.class);

    private ServicesProvider() {
        // hidden constructor
    }

    public static final WebTarget getWebTarget(final String servicesProvider) {
        Preconditions.checkNotNull(servicesProvider, "Services provider URL");
        final StringBuilder servicesStringBuilder = new StringBuilder(servicesProvider);
        if ((!servicesProvider.endsWith(PUTWSS)) && (!servicesProvider.endsWith(PUTWSS + "/"))) {
            if (!servicesProvider.endsWith("/")) {
                servicesStringBuilder.append("/");
            }
            servicesStringBuilder.append(PUTWSS);
        }
        final Client client = createClient();
        if (client != null) {
            return client.target(UriBuilder.fromUri(servicesStringBuilder.toString()).build());
        }
        throw new IllegalStateException("Service provider cannot be initialized!");
    }

    private static Client createClient() {
        try {
            final ClientConfig config = new ClientConfig();
            final SSLContext ctx = SSLContext.getInstance("TLS");
            ctx.init(null, new TrustManager[] { getDefaultTrustManager() }, new SecureRandom());
            return ClientBuilder.newBuilder().withConfig(config).sslContext(ctx).build();
        } catch (Exception e) {
            LOGGER.error(e.getMessage(), e);
        }
        return null;
    }

    private static X509TrustManager getDefaultTrustManager() {
        try {
            final TrustManagerFactory factory = TrustManagerFactory.getInstance(TrustManagerFactory
                    .getDefaultAlgorithm());
            factory.init((KeyStore) null);
            return Iterables
                    .getFirst(Iterables.filter(Arrays.asList(factory.getTrustManagers()),
                            X509TrustManager.class), null);
        } catch (NoSuchAlgorithmException | KeyStoreException e) {
            LOGGER.error(e.getMessage(), e);
        }
        return null;
    }
}
