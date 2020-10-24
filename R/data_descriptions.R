#' Full data from Blanford et al (2012)
#'
#' Complete data from the publication: Blanford S, Jenkins NE, Read AF, Thomas
#' MB (2012) Evaluating the lethal and pre-lethal effects of a range of fungi
#' against adult Anopheles stephensi mosquitoes. Malaria Journal. 11:365
#' \href{https://doi.org/10.1186/1475-2875-11-365}{doi}
#'
#' \describe{
#'   \item{block}{experimental block within experiment (1 - 5)}
#'   \item{treatment}{experimental treatment}
#'   \item{replicate cage}{replicate cage within treatment (1 - 4)}
#'   \item{day}{time post-infection (days)}
#'   \item{censor}{'1' censored, '0' died }
#'   \item{d}{an indicator variable; '0' censored, '1' died}
#'   \item{inf}{'0' uninfected treatement, '1' infected treatment}
#'   \item{t}{time post-infection (days)}
#'   \item{fq}{frequency of individuals}
#'   }
#'
#' @source Simon Blanford and Matthew Thomas generously provided and allowed the
#'   release of these data
#' @examples
#' head(data_blanford)
"data_blanford"


#' Full data from Parker et al (2014)
#'
#' Complete data on the survival of adult female aphids exposed or
#' unexposed to fungal infection.
#'
#' @return A dataframe
#' \describe{
#' \item{Genotype}{names of host genotypes}
#' \item{SD}{Spore Dose, concentration of fungal spores hosts were exposed to
#' (spores/mm^2)}
#' \item{Fecundity}{total number of offspring produced by hosts over lifetime}
#' \item{Sporulation}{'1' visual signs of sporulation at host death,
#' '0' no signs of sporulation}
#' \item{Status}{'0' censored, '1' died}
#' \item{Time}{time of death (days)}
#' \item{dose}{dose treatments on ordinal scale of 1-3, controls 0}
#' \item{censored}{'1' censored, '0' died}
#' \item{d}{death indicator: '1' died, '0' censored}
#' \item{t}{time of death (days)}
#' \item{g}{infection treatment indicator; '1' infected, '0' uninfected}
#' }
#' @source Parker BJ, Garcia JR, Gerardo NM (2014)
#' Genetic variation in resistance and fecundity tolerance in a natural
#' host-pathogen interaction.
#' Evolution 68: 2421-2429 \href{https://doi.org/10.1111/evo.12418}{doi}
#'
#' The full dataset is available at Dryad
#' \url{https://doi.org/10.5061/dryad.24gq7}
#' @examples
#' head(data_parker)
"data_parker"


#' A subset of data from Lorenz & Koella (2011)
#'
#' Data on the longevity of 256 adult female mosquitoes and the number of
#' pathogen spores they harboured at the time of their death.
#'
#' These are the Lorenz & Koella data analysed in Agnew (2019)
#' \href{https://doi.org/10.1101/530709}{doi}
#' @return A dataframe
#' \describe{
#' \item{Infectious.dose}{Number of spores larvae
#' were exposed to (spores/larva)}
#' \item{Food}{food treatment: '50' or '100'}
#' \item{Sex}{sex of mosquito: 'F' female}
#' \item{Spore.Count}{number of spores harboured by mosquito at time of death}
#' \item{t}{time of death to nearest half day}
#' \item{censored}{'1' censored, '0' died}
#' \item{d}{death indicator: '1' died during experiment, '0' right-censored}
#' \item{g}{infection treatment indicator: '0' uninfected, '1' infected}
#' }
#' @source Lorenz LM & Koella JC (2011) The microsporidian parasite
#' Vavraia culicis as a potential late life-acting control agent of malaria.
#' Evolutionary
#' Applications 4: 783-790
#' \href{https://doi.org/10.1111/j.1752-4571.2011.00199.x}{doi}
#'
#' The full dataset is available at Dryad
#' \url{https://doi.org/10.5061/dryad.2s231}
#' @examples
#' head(data_lorenz)
"data_lorenz"


#' Simulated recovery data
#'
#' Simulated data allowing for recovery from infection,
#' where recovered individuals experience the same background mortality
#' as uninfected individuals.
#' @return A dataframe
#' \describe{
#' \item{columns 1-2}{indicator variables identifying individuals
#' in the control treatment that died (control.d = 1)
#' or were censored (control.c = 1)}
#' \item{columns 3-4}{indicator variables identifying individuals
#' in the infected treatment that died (infected.d = 1)
#' or were censored (infected.c = 1) while still infected}
#' \item{columns 5-6}{indicator variables indenifying individuals
#' in the infected treatment that died (recovered.d = 1)
#' or were censored (recovered.c = 1) after recovering from infection}
#' \item{censor}{'0' not censored, '1' censored}
#' \item{d}{death indicator: '1' died, '0'}
#' \item{t}{time of the event (death or censoring)}
#' \item{fq}{frequence or number of individuals associated with the event
#' of death or censoring at time t}
#' }
#' @examples
#' head(recovery_data) ; tail(recovery_data)
"recovery_data"


#' Simulated recovery data, with no background mortality
#'
#' Simulated data allowing for recovery from infection,
#' when there is no background mortality.
#' @return A dataframe
#' \describe{
#' \item{columns 1-2}{indicator variables identifying individuals
#' in the control treatment that died (control.d = 1)
#' or were censored (control.c = 1)}
#' \item{columns 3-4}{indicator variables identifying individuals
#' in the infected treatment that died (infected.d = 1)
#' or were censored (infected.c = 1) while still infected}
#' \item{columns 5-6}{indicator variables indenifying individuals
#' in the infected treatment that died (recovered.d = 1)
#' or were censored (recovered.c = 1) after recovering from infection}
#' \item{censor}{'0' not censored, '1' censored}
#' \item{d}{death indicator: '1' died, '0'}
#' \item{t}{time of the event (death or censoring)}
#' \item{fq}{frequence or number of individuals associated with the event
#' of death or censoring at time t}
#' }
#' @examples
#' head(recovery_data_II) ; tail(recovery_data_II)
"recovery_data_II"




