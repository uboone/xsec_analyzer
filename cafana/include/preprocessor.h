/**
 * @file preprocessor.h
 * @brief Header file for preprocessor macros that streamline variable
 * declarations.
 * @author mueller@fnal.gov
*/
#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H
#include "sbnana/CAFAna/Core/MultiVar.h"

/**
 * @brief Preprocessor wrapper for looping over reco interactions and
 * broadcasting a function over the reco interactions.
 * @details This macro declares a lambda function that broadcasts a function
 * over the reco interactions and returns a vector of the results of the
 * function.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each reco interaction
 * passing the cut SEL.
 */
#define SPINEVAR_RR(VAR,SEL)                             \
    [](const caf::SRSpillProxy* sr)                      \
    {                                                    \
        std::vector<double> var;                         \
        bool is_mc(sr->ndlp_true != 0);                  \
        for(auto const& i : sr->dlp)                     \
        {                                                \
            if(SEL(i) && (i.match.size() > 0 || !is_mc)) \
                var.push_back(VAR(i));                   \
        }                                                \
        return var;                                      \
    }

/**
 * @brief Preprocessor wrapper for looping over reco interactions and
 * broadcasting a function over the matching true interactions.
 * @details This macro declares a lambda function that broadcasts a function
 * over the reco interactions and returns a vector of the results of the
 * function for each true interaction.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each true interaction that
 * is matched to by a reco interaction passing the cut SEL.
 */
#define SPINEVAR_RT(VAR,SEL)                                                              \
    [](const caf::SRSpillProxy* sr)                                                       \
    {                                                                                     \
        std::vector<double> var;                                                          \
        bool is_mc(sr->ndlp_true != 0);                                                   \
        for(auto const& i : sr->dlp)                                                      \
        {                                                                                 \
            if(SEL(i) && (i.match.size() > 0 || !is_mc))                                  \
                var.push_back(i.match.size() > 0 ? VAR(sr->dlp_true[i.match[0]]) : -1.0); \
        }                                                                                 \
        return var;                                                                       \
 }

/**
 * @brief Preprocessor wrapper for looping over true interactions and
 * broadcasting a function over the matching reco interactions.
 * @details This macro declares a lambda function that broadcasts a function
 * over the true interactions and returns a vector of the results of the
 * function for each reco interaction.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each reco interaction that
 * is matched to by a true interaction passing the cut SEL.
 */
#define SPINEVAR_TR(VAR,SEL)                         \
    [](const caf::SRSpillProxy* sr)                  \
    {							                     \
        std::vector<double> var;			         \
        for(auto const& i : sr->dlp_true)	         \
        {                                            \
            if(SEL(i) && i.match.size() > 0)         \
            var.push_back(VAR(sr->dlp[i.match[0]])); \
        }                                            \
        return var;                                  \
    }

/**
 * @brief Preprocessor wrapper for looping over true interactions and
 * broadcasting a function over the true interactions.
 * @details This macro declares a lambda function that broadcasts a function
 * over the true interactions and returns a vector of the results of the
 * function.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each true interaction
 * passing the cut SEL.
 */
#define SPINEVAR_TT(VAR,SEL)                 \
    [](const caf::SRSpillProxy* sr)          \
    {                                        \
        std::vector<double> var;             \
        for(auto const& i : sr->dlp_true)    \
        {                                    \
            if(SEL(i) && i.match.size() > 0) \
                var.push_back(VAR(i));       \
        }                                    \
        return var;                          \
    }

/**
 * @brief Preprocessor wrapper for looping over slices and broadcasting
 * a function over the slices.
 * @details This macro declares a lambda function that broadcasts a function
 * over the slices and returns a vector of the results of the function.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each true interaction
 * passing the cut SEL.
 */

#define PANDORAVAR(VAR,SEL)                  \
    [](const caf::SRSpillProxy* sr)          \
    {                                        \
        std::vector<double> var;             \
        for(auto const& i : sr->slc)         \
        {                                    \
            if(SEL(i))                       \
                var.push_back(VAR(i));       \
        }                                    \
        return var;                          \
    }

#define PANDORAMULTIVAR(VAR,SEL)              \
    [](const caf::SRSpillProxy* sr)          \
    {                                        \
        std::vector<double> var;       \
        for(auto const& i : sr->slc)         \
        {                                 \
            auto res = VAR(i);                                  \
            if(SEL(i))                      \
                var.insert(var.end(), res.begin(), res.end());                       \
        }                                    \
        return var;                          \
    }
#endif // PREPROCESSOR_H