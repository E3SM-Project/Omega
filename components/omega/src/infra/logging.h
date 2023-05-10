// Omega logging header file.

#ifndef OMEGA_LOG_H
#define OMEGA_LOG_H

#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/null_sink.h>

auto CONSOLE = spdlog::stdout_logger_mt("*");

#if defined(OMEGA_LOG_OUTFILE)
    auto LOGFILE = spdlog::basic_logger_mt("F", "OMEGA_LOG_OUTFILE");
#else
    auto LOGFILE = spdlog::null_logger_mt("F");
#endif

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
#    define LOGGER_TRACE(logger, ...) SPDLOG_LOGGER_TRACE(logger, __VA_ARGS__)
#    define LOG_TRACE(...) LOGGER_TRACE(CONSOLE, __VA_ARGS__)
#else
#    define LOGGER_TRACE(logger, ...) (void)0
#    define LOG_TRACE(...) (void)0
#endif

#if OMEGA_LOG_ACTIVE_LEVEL <= OMEGA_LOG_LEVEL_DEBUG
#    define LOGGER_DEBUG(logger, ...) SPDLOG_LOGGER_DEBUG(logger, __VA_ARGS__)
#    define LOG_DEBUG(...) LOGGER_DEBUG(CONSOLE, __VA_ARGS__)
#else
#    define LOGGER_DEBUG(logger, ...) (void)0
#    define LOG_DEBUG(...) (void)0
#endif

#if OMEGA_LOG_ACTIVE_LEVEL <= OMEGA_LOG_LEVEL_INFO
#    define LOGGER_INFO(logger, ...) SPDLOG_LOGGER_INFO(logger, __VA_ARGS__)
#    define LOG_INFO(...) LOGGER_INFO(CONSOLE, __VA_ARGS__)
#else
#    define LOGGER_INFO(logger, ...) (void)0
#    define LOG_INFO(...) (void)0
#endif

#if OMEGA_LOG_ACTIVE_LEVEL <= OMEGA_LOG_LEVEL_WARN
#    define LOGGER_WARN(logger, ...) SPDLOG_LOGGER_WARN(logger, __VA_ARGS__)
#    define LOG_WARN(...) LOGGER_WARN(CONSOLE, __VA_ARGS__)
#else
#    define LOGGER_WARN(logger, ...) (void)0
#    define LOG_WARN(...) (void)0
#endif

#if OMEGA_LOG_ACTIVE_LEVEL <= OMEGA_LOG_LEVEL_ERROR
#    define LOGGER_ERROR(logger, ...) SPDLOG_LOGGER_ERROR(logger, __VA_ARGS__)
#    define LOG_ERROR(...) LOGGER_ERROR(CONSOLE, __VA_ARGS__)
#else
#    define LOGGER_ERROR(logger, ...) (void)0
#    define LOG_ERROR(...) (void)0
#endif

#if OMEGA_LOG_ACTIVE_LEVEL <= OMEGA_LOG_LEVEL_CRITICAL
#    define LOGGER_CRITICAL(logger, ...) SPDLOG_LOGGER_CRITICAL(logger, __VA_ARGS__)
#    define LOG_CRITICAL(...) LOGGER_CRITICAL(CONSOLE, __VA_ARGS__)
#else
#    define LOGGER_CRITICAL(logger, ...) (void)0
#    define LOG_CRITICAL(...) (void)0
#endif

#endif // OMEGA_LOG_H
