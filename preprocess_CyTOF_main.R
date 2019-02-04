# A more automatic version of the script for pre-processing of CyTOF data - but without GUI usage as this is too slow

#### NOTES - please read ####
# Because of the inherent need to iterate and plot data in the pre-processing of CyTOF-data, this script will have to be run 
# in different modes - the first part may be run 'once' in the beginning, but the rest requires you to input which batch to work with
# and then you can process these one at a time.



#### Sourcing necessary functions
source('pregating_standalone.R')
source('normalization_standalone.R')

# Set wd
work.dir = '/work/directory/'

# First, specify files to work on
batch1_files <- list.files(file.path(work.dir, 'Raw_data', 'batch1'), pattern = '\\.FCS', full.names = T)
batch2_files <- list.files(file.path(work.dir, 'Raw_data', 'batch2'), pattern = '\\.FCS', full.names = T)

# Specify whether to concat in the order given in the time variable for the FCS files - this is T by default (T), but sometimes you
# have batches run over multiple days and it does not consider dates...
# In those cases, make sure your batchX_files contains the individual FCS files in the order, you wish to concatenate in.
concat_by_time <- c(T, T)
n_batches = length(grep('^batch[0-9]+_files$', ls()))

# Quick check for consistency
if (length(concat_by_time) != n_batches) {
  stop('Error: Your number of batches does not match the length of the concat_by_time variable!\n')
}




#### Processing for ALL files - give this a lot of memory! For 8 batches I use a whole node with 350 gb memory - it's fast however - it only takes about three hours ####

#### Step 1: Concatenation of files ####
concat_batches <- list()
for (batch in 1:n_batches) {

  # Check if batch exists
  batch_name = paste0('batch', batch)
  if (exists(paste0(batch_name, '_files'))) {
    concat_batches[[batch_name]] <-  CATALYST::concatFCS(eval(parse(text=paste0(batch_name, '_files'))), by_time = concat_by_time[batch])
  } else {
    stop('Error: Your batch names are not correct. They should be batch1_files, batch2_files, ..., batchN_files.\n')
  }

}



#### Step 2: Normalization of expression data using MATLAB normalizer ####
# This requires us to see a plot and specify boundaries - first the plots - now look at these and decide on the bondaries you want relative to the reference set
if (!dir.exists(file.path(work.dir, 'normalization_plots'))) {dir.create(file.path(work.dir, 'normalization_plots'))}
#lapply(concat_batches, plot_beads, xleft = 4, xright = 9, ybottom = -0.2, ytop = 4, save = T, wd = file.path(work.dir, 'normalization_plots'))


# After looking at the plots, set the boundaries in the below command and go - for some reason plotting will not function here...
normed_batches <- normalizeFCS_cbp(concat_batches, bead_type = 'Fluidigm', wd = file.path(work.dir, 'normalization_plots'),
                                   beads_gates = list(list("Ce140Di" = list(x = c(4, 9), y = c(-0.5, 4)),
                                                           "Eu151Di" = list(x = c(4, 9), y = c(-0.5, 4)),
                                                           "Eu153Di" = list(x = c(4, 9), y = c(-0.5, 4)),
                                                           "Ho165Di" = list(x = c(4, 9), y = c(-0.5, 4)),
                                                           "Lu175Di" = list(x = c(4.2, 9), y = c(-0.5, 3.5))),

                                                      list("Ce140Di" = list(x = c(4, 9), y = c(-0.5, 4.2)),
                                                           "Eu151Di" = list(x = c(4, 9), y = c(-0.5, 4.2)),
                                                           "Eu153Di" = list(x = c(4, 9), y = c(-0.5, 4.2)),
                                                           "Ho165Di" = list(x = c(4, 9), y = c(-0.5, 4.2)),
                                                           "Lu175Di" = list(x = c(4, 9), y = c(-0.5, 4.2)))))


# Since the following steps are conducted one at a time for the fcs files, a good solution to avoid memory issues is to
# save the normed data as multiple R objects, clear R's memory and read the normed data in, one file at a time
# Here, I save all the normed objects one by one, delete the whole set and load them in one-by-one for downstream processing
if (!dir.exists(file.path(work.dir, 'normalized'))) {dir.create(file.path(work.dir, 'normalized'))}

for (batch in 1:n_batches) {
  saveRDS(normed_batches[[batch]], file = file.path(work.dir, 'normalized', paste0('normalized_batch', batch, '.rds')))
}

# Remove whole set to clear memory
rm(concat_batches, normed_batches)




#### Step 3: Calculation of spillover matrix - from CATALYST ####
# Read single stained beads file
ss_beads <- flowCore::read.FCS(file.path(work.dir, 'single_stained_beads.fcs'), truncate_max_range = F, normalization = F)

# Specify the mass channels stained for
bc_ms <- c(89, 114:115, 141:156, 158:176, 198, 209)

# Debarcode and compute spill matrix
re <- CATALYST::assignPrelim(x = ss_beads, y = bc_ms)
re <- CATALYST::estCutoffs(x = re)
re <- CATALYST::applyCutoffs(x = re)
spill_mat <- CATALYST::computeSpillmat(x = re, method = 'default', interactions = 'default')

# Check the plot and see if you have a nice estimations
#plotSpillmat(bc_ms=c(89, 114:115, 141:156, 158:176, 198, 209), SM=spill_mat, plotly=FALSE, isotope_list = c(CATALYST::isotope_list, list(BCKG=190)))

# Save the spillover matrix
if (!dir.exists(file.path(work.dir, 'compensation'))) {dir.create(file.path(work.dir, 'compensation'))}
saveRDS(spill_mat, file = file.path(work.dir, 'compensation', paste0('spillover_matrix.rds')))


# Cleaning memory
rm(ss_beads, bc_ms, re, spill_mat)




#### Processing for individual files ####
# You can out-comment all the rest if you wish to run the above at first in one go.

#### Step 1: Removal of doublets, beads, dead cells, debris etc. ####

# This process also requires some plotting to see what is going on - first choose the channels to use in gating
channels_to_plot <- list(c('Ce140Di', 'Ir193Di'), c('Ir191Di', 'Ir193Di'), c('Event_length', 'Ir193Di'), c('Pt194Di', 'Pt195Di'))

# Set the cutoffs you wish to use in the four pregating steps - for ease you set these for all files in the same script
pre_gates <- list(list(list(x = c(-0.2, 3.8), y = c(5.2, 9.5)),
                       list(x = c(5.3, 6.3), y = c(5.9, 6.8)),
                       list(x = c(10, 30), y = c(5.9, 6.8)),
                       list(x = c(-0.2, 3.8), y = c(-0.2, 2.1))),

                  list(list(x = c(-0.2, 3.8), y = c(5.2, 9.5)),
                       list(x = c(5.25, 6.2), y = c(5.9, 6.8)),
                       list(x = c(10, 30), y = c(5.9, 6.8)),
                       list(x = c(-0.2, 3.8), y = c(-0.2, 2.5))))


# Check that you set enough gates
if (length(pre_gates) != n_batches) {
  stop('Error: Your number of batches does not match the length of the pre_gates list!\n')
}


# Plot strategy - this is meant to be an iterative process in which you plot until you're satisfied with the settings
# Because this requires loading etc. it's hard to make a good loop for it. As a consequence, I have just made
# it possible to set a batch numnber to work on...
working_batch = 1




# Now we work with that one batch
# Loading batch data
normed_batch <- readRDS(file.path(work.dir, 'normalized', paste0('normalized_batch', working_batch, '.rds')))

# Check gates OK format
if (length(pre_gates[[working_batch]]) != length(channels_to_plot)) {
  stop('Error: Your set gates do not have the right format to fit with the specified plotting channels.\n')
}

# Visualize strategy
if (!dir.exists(file.path(work.dir, 'pregating_plots'))) {dir.create(file.path(work.dir, 'pregating_plots'))}
#visualize_strategy(normed_batch, paste0('batch', working_batch), wd = file.path(work.dir, 'pregating_plots'), channels_to_plot, pre_gates[[working_batch]])


# Apply strategy
pregated_batch <- apply_strategy(normed_batch, paste0('batch', working_batch), wd = file.path(work.dir, 'pregating_plots'), channels_to_plot, pre_gates[[working_batch]])
if (!dir.exists(file.path(work.dir, 'pregated_data'))) {dir.create(file.path(work.dir, 'pregated_data'))}
saveRDS(pregated_batch, file = file.path(work.dir, 'pregated_data', paste0('pregated_batch', working_batch, '.rds')))

# Clean up
rm(normed_batch)



#### Step 2: Compensation of spillover - from CATALYST ####
# Read spillover matrix
spill_mat <- readRDS(file.path(work.dir, 'compensation', paste0('spillover_matrix.rds')))

# Actual compensation step
comped_data <- CATALYST::compCytof(x = pregated_batch, y = spill_mat, method="nnls", isotope_list = c(CATALYST::isotope_list, list(BCKG=190)))
# Takes about 1,5 hours on 10 cores with 100 gb memory PER file...

# We clean up and save the comped_data object
rm(pregated_batch, spill_mat)
saveRDS(comped_data, file = file.path(work.dir, 'compensation', paste0('compensated_batch', working_batch, '.rds')))




#### Step 3: Debarcoding ####
# Specify debarcoding scheme
barcodes <- lapply(1:n_batches, function(x) {openxlsx::read.xlsx(file.path(work.dir, 'debarcoding_scheme.xlsx'), sheet = paste0('batch', x), rowNames = T)})


# Set the cutoffs you wish to use in the four pregating steps - for ease you set these for all files in the same script
sep_cutoffs_batches <- list(c(0.32,0.32,0.32,0.32,0.33,
                              0.30,0.29,0.28,0.27,0.28,
                              0.30,0.28,0.30,0.32,0.30,
                              0.29,0.35,0.34,0.32,0.29),

                            c(0.27,0.33,0.25,0.25,0.29,
                              0.29,0.23,0.25,0.25,0.22,
                              0.27,0.23,0.25,0.23,0.27,
                              0.27,0.23,0.29,0.26,0.25))


# Assigning each event a preliminary barcode ID
re0 <- CATALYST::assignPrelim(x = comped_data, y = barcodes[[working_batch]], verbose=FALSE)

# Estimate a sample-specific cutoff - reviewing the yield plots is still advisable
re1 <- CATALYST::estCutoffs(x=re0)

# Use population-specific cutoffs
re <- CATALYST::applyCutoffs(x = re1, sep_cutoffs = sep_cutoffs_batches[[working_batch]])

# Yield plot
if (!dir.exists(file.path(work.dir, 'debarcoding'))) {dir.create(file.path(work.dir, 'debarcoding'))}
pdf(file.path(work.dir, 'debarcoding', paste0('batch', working_batch, '_yields.pdf')), height = 5, width = 10)
CATALYST::plotYields(x=re, which = 0, plotly = FALSE)
dev.off()

# Make plots of utoffs for all samples
pdf(file.path(work.dir, 'debarcoding', paste0('event_plots_batch', working_batch, '.pdf')), height = 5, width = 10)
CATALYST::plotEvents(x = re, which = 'all', n_events = 500)
dev.off()

# Save objects to files
saveRDS(re0, file = file.path(work.dir, 'debarcoding', paste0('debarcoding_prelim_batch', working_batch, '_1.rds')))
saveRDS(re1, file = file.path(work.dir, 'debarcoding', paste0('debarcoding_prelim_batch', working_batch, '_2.rds')))




# Write to files
if (!dir.exists(file.path(work.dir, 'debarcoding', paste0('batch', working_batch)))) {dir.create(file.path(work.dir, 'debarcoding', paste0('batch', working_batch)))}
CATALYST::outFCS(x = re, y = comped_data, out_path = file.path(work.dir, 'debarcoding', paste0('batch', working_batch)))


# Make excel report
debarcoding_report <- cbind(unname(table(re@bc_ids)), unname(table(re@bc_ids))/length(re@bc_ids), unname(re@sep_cutoffs[names(table(re@bc_ids))]))
colnames(debarcoding_report) <- c('Cell count', 'Percentage', 'Debarcoding cut-off'); rownames(debarcoding_report) <- names(table(re@bc_ids))
openxlsx::write.xlsx(debarcoding_report, file = file.path(work.dir, 'debarcoding', paste0('debarcoding_report_batch', working_batch, '.xlsx')), row.names=TRUE, sheetName = paste0('batch', working_batch))
