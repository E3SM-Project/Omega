// Omega logging header file.

#ifndef OMEGA_LOG_H
#define OMEGA_LOG_H

#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/null_sink.h>

#define OMEGA_LOG_LEVEL_TRACE SPDLOG_LEVEL_TRACE
#define OMEGA_LOG_LEVEL_DEBUG SPDLOG_LEVEL_DEBUG
#define OMEGA_LOG_LEVEL_INFO SPDLOG_LEVEL_INFO
#define OMEGA_LOG_LEVEL_WARN SPDLOG_LEVEL_WARN
#define OMEGA_LOG_LEVEL_ERROR SPDLOG_LEVEL_ERROR
#define OMEGA_LOG_LEVEL_CRITICAL SPDLOG_LEVEL_CRITICAL
#define OMEGA_LOG_LEVEL_OFF SPDLOG_LEVEL_OFF

#if !defined(OMEGA_LOG_ACTIVE_LEVEL)
#    define OMEGA_LOG_ACTIVE_LEVEL OMEGA_LOG_LEVEL_INFO
#endif

#if OMEGA_LOG_ACTIVE_LEVEL <= OMEGA_LOG_LEVEL_TRACE
#    define OMEGA_LOGGER_TRACE(logger, ...) SPDLOG_LOGGER_TRACE(logger, __VA_ARGS__)
#    define OMEGA_LOG_TRACE(...) OMEGA_LOGGER_TRACE(spdlog::default_logger(), __VA_ARGS__)
#else
#    define OMEGA_LOGGER_TRACE(logger, ...) (void)0
#    define OMEGA_LOG_TRACE(...) (void)0
#endif

#if OMEGA_LOG_ACTIVE_LEVEL <= OMEGA_LOG_LEVEL_DEBUG
#    define OMEGA_LOGGER_DEBUG(logger, ...) SPDLOG_LOGGER_DEBUG(logger, __VA_ARGS__)
#    define OMEGA_LOG_DEBUG(...) OMEGA_LOGGER_DEBUG(spdlog::default_logger(), __VA_ARGS__)
#else
#    define OMEGA_LOGGER_DEBUG(logger, ...) (void)0
#    define OMEGA_LOG_DEBUG(...) (void)0
#endif

#if OMEGA_LOG_ACTIVE_LEVEL <= OMEGA_LOG_LEVEL_INFO
#    define OMEGA_LOGGER_INFO(logger, ...) SPDLOG_LOGGER_INFO(logger, __VA_ARGS__)
#    define OMEGA_LOG_INFO(...) OMEGA_LOGGER_INFO(spdlog::default_logger(), __VA_ARGS__)
#else
#    define OMEGA_LOGGER_INFO(logger, ...) (void)0
#    define OMEGA_LOG_INFO(...) (void)0
#endif

#if OMEGA_LOG_ACTIVE_LEVEL <= OMEGA_LOG_LEVEL_WARN
#    define OMEGA_LOGGER_WARN(logger, ...) SPDLOG_LOGGER_WARN(logger, __VA_ARGS__)
#    define OMEGA_LOG_WARN(...) OMEGA_LOGGER_WARN(spdlog::default_logger(), __VA_ARGS__)
#else
#    define OMEGA_LOGGER_WARN(logger, ...) (void)0
#    define OMEGA_LOG_WARN(...) (void)0
#endif

#if OMEGA_LOG_ACTIVE_LEVEL <= OMEGA_LOG_LEVEL_ERROR
#    define OMEGA_LOGGER_ERROR(logger, ...) SPDLOG_LOGGER_ERROR(logger, __VA_ARGS__)
#    define OMEGA_LOG_ERROR(...) OMEGA_LOGGER_ERROR(spdlog::default_logger(), __VA_ARGS__)
#else
#    define OMEGA_LOGGER_ERROR(logger, ...) (void)0
#    define OMEGA_LOG_ERROR(...) (void)0
#endif

#if OMEGA_LOG_ACTIVE_LEVEL <= OMEGA_LOG_LEVEL_CRITICAL
#    define OMEGA_LOGGER_CRITICAL(logger, ...) SPDLOG_LOGGER_CRITICAL(logger, __VA_ARGS__)
#    define OMEGA_LOG_CRITICAL(...) OMEGA_LOGGER_CRITICAL(spdlog::default_logger(), __VA_ARGS__)
#else
#    define OMEGA_LOGGER_CRITICAL(logger, ...) (void)0
#    define OMEGA_LOG_CRITICAL(...) (void)0
#endif

#include "YAKL.h"

typedef yakl::Array<double,2,yakl::memDevice> doub2d;

template<>
struct fmt::formatter<doub2d> : fmt::formatter<std::string>
{
    auto format(doub2d my, format_context &ctx) -> decltype(ctx.out())
    {
        return fmt::format_to(ctx.out(), "[doub2d ={}]", my.label());
    }
};

#endif // OMEGA_LOG_H
