#' dataset of Long et al. (2020)
#'
#'
#' @format A `data.frame` for gamma distributed outcome with 5 studies:
#' \describe{
#'  \item{Study}{Study ID}
#'  \item{mi}{Mean of intensive care unit (ICU) length of stay (LoS)}
#'  \item{sdi}{Standard deviation of ICU LoS}
#'  \item{ni}{Number of subjects}
#'  \item{trt}{1 is Treatment, an 0 is control}
#' }
#'
#' @source Long, R., Tian, J., Wu, S., Li, Y., Yang, X., & Fei, J. (2020). Clinical efficacy of surgical versus conservative treatment for multiple rib fractures: a meta-analysis of randomized controlled trials. *International Journal of Surgery*, 83, 79-88.
#'   Data used under the terms of the Creative Commons Attribution 4.0 International License (CC BY 4.0):
#'   <https://creativecommons.org/licenses/by/4.0/>
"long2020"


#' dataset of Chu et al. (2020)
#'
#'
#' @format A `data.frame` for binary distributed outcome with 7 studies:
#' \describe{
#'  \item{Study}{Study ID}
#'  \item{x1}{Events of Covid-19 for group of further distnce}
#'  \item{n1}{Number of subjects of Covid-19 for group of further distnce}
#'  \item{x0}{Events of Covid-19 for group of shorter distnce}
#'  \item{n0}{Number of subjects of Covid-19 for group of shorter distnce}
#' }
#'
#' @source Chu, D. K., Akl, E. A., Duda, S., Solo, K., Yaacoub, S., Schünemann, H. J., ... & Reinap, M. (2020). Physical distancing, face masks, and eye protection to prevent person-to-person transmission of SARS-CoV-2 and COVID-19: a systematic review and meta-analysis. *The lancet*, 395(10242), 1973-1987.
#'   Data used under the terms of the Creative Commons Attribution IGO 3.0 License (CC BY 3.0 IGO):
#'   <https://creativecommons.org/licenses/by/3.0/igo/>
"chu2020"


#' dataset of Rutter et al. (2021)
#'
#'
#' @format A `data.frame` for Poisson distributed outcome with 11 studies:
#' \describe{
#'  \item{Study}{Study ID}
#'  \item{Year}{Publishing year}
#'  \item{Sex}{Sex (Female or Male)}
#'  \item{Events}{Number of events}
#'  \item{Time}{Total observation time}
#' }
#'
#' @source Rutter, M., Bowley, J., Lanyon, P. C., Grainge, M. J., & Pearce, F. A. (2021). A systematic review and meta-analysis of the incidence rate of Takayasu arteritis. *Rheumatology*, 60(11), 4982-4990.
#'   Data used under the terms of the Creative Commons Attribution 4.0 International License (CC BY 4.0):
#'   <https://creativecommons.org/licenses/by/4.0/>
"rutter2021"


