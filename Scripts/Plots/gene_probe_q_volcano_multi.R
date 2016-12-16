# Scatterplot of median RMA-normalized values
q_data <- data.frame();
for (GSE in dbs) {
    q_data <- rbind.data.frame(
        q_data,
        cbind.data.frame(
            "Slope" = regr_values_rma[[GSE]][,"coeff"],
            "Q" = -log10(regr_values_rma[[GSE]][,"q"]),
            "Dataset" = GSE
        )
    );
}
png("GPL570_rma_q.png", width = 3000, height = 1800, res = 300);
p1 <- (
    ggplot(
        data = q_data,
        mapping = aes(x = Slope, y = Q, colour = Dataset))
    + geom_point(alpha = 0.2)
    + geom_hline(yintercept = 2, linetype = 2)
    + ggtitle("Significant probe expression")
    + xlab("Linear regression slope")
    + ylab("-log10(q) (FDR adjusted)")
    + guides(colour = guide_legend(override.aes = list(alpha = 1)))
    + xlim(c(-1e-6, 1e-6))
    # + ylim(c(4,6))
    # + scale_x_datetime(date_labels = "%Y-%m-%d", date_breaks = "4 month")
    # + theme(
    #     axis.text.x = element_text(angle = 90),
    #     panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
    #     panel.grid.minor = element_blank(), 
    #     panel.grid.major = element_blank(),
    #     plot.background = element_rect(fill = "transparent",colour = NA)
    # )
);
print(p1);
dev.off();