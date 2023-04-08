//
// Copyright(c) 2015 Gabi Melman.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

// spdlog usage example

#include <cstdio>
#include <chrono>

void load_levels_example();
void stdout_logger_example();
void basic_example();
void rotating_example();
void daily_example();
void callback_example();
void async_example();
void binary_example();
void vector_example();
void stopwatch_example();
void trace_example();
void multi_sink_example();
void user_defined_example();
void err_handler_example();
void syslog_example();
void udp_example();
void custom_flags_example();
void file_events_example();
void replace_default_logger_example();

#include "spdlog/spdlog.h"
#include "spdlog/cfg/env.h"  // support for loading levels from the environment variable
#include "spdlog/fmt/ostr.h" // support for user defined types

int main(int, char *[])
{
    // Log levels can be loaded from argv/env using "SPDLOG_LEVEL"
    load_levels_example();

    spdlog::info("Welcome to spdlog version {}.{}.{}  !", SPDLOG_VER_MAJOR, SPDLOG_VER_MINOR, SPDLOG_VER_PATCH);

    spdlog::warn("Easy padding in numbers like {:08d}", 12);
    spdlog::critical("Support for int: {0:d};  hex: {0:x};  oct: {0:o}; bin: {0:b}", 42);
    spdlog::info("Support for floats {:03.2f}", 1.23456);
    spdlog::info("Positional args are {1} {0}..", "too", "supported");
    spdlog::info("{:>8} aligned, {:<8} aligned", "right", "left");

    // Runtime log levels
    spdlog::set_level(spdlog::level::info); // Set global log level to info
    spdlog::debug("This message should not be displayed!");
    spdlog::set_level(spdlog::level::trace); // Set specific logger's log level
    spdlog::debug("This message should be displayed..");

    // Customize msg format for all loggers
    spdlog::set_pattern("[%H:%M:%S %z] [%^%L%$] [thread %t] %v");
    spdlog::info("This an info message with custom format");
    spdlog::set_pattern("%+"); // back to default format
    spdlog::set_level(spdlog::level::info);

    // Backtrace support
    // Loggers can store in a ring buffer all messages (including debug/trace) for later inspection.
    // When needed, call dump_backtrace() to see what happened:
    spdlog::enable_backtrace(10); // create ring buffer with capacity of 10  messages
    for (int i = 0; i < 100; i++)
    {
        spdlog::debug("Backtrace message {}", i); // not logged..
    }
    // e.g. if some error happened:
    spdlog::dump_backtrace(); // log them now!

    try
    {
        stdout_logger_example();
        basic_example();
        rotating_example();
        daily_example();
        callback_example();
        async_example();
        binary_example();
        vector_example();
        multi_sink_example();
        user_defined_example();
        err_handler_example();
        trace_example();
        stopwatch_example();
        udp_example();
        custom_flags_example();
        file_events_example();
        replace_default_logger_example();

        // Flush all *registered* loggers using a worker thread every 3 seconds.
        // note: registered loggers *must* be thread safe for this to work correctly!
        spdlog::flush_every(std::chrono::seconds(3));

        // Apply some function on all registered loggers
        spdlog::apply_all([&](std::shared_ptr<spdlog::logger> l) { l->info("End of example."); });

        // Release all spdlog resources, and drop all loggers in the registry.
        // This is optional (only mandatory if using windows + async log).
        spdlog::shutdown();
    }

    // Exceptions will only be thrown upon failed logger or sink construction (not during logging).
    catch (const spdlog::spdlog_ex &ex)
    {
        std::printf("Log initialization failed: %s\n", ex.what());
        return 1;
    }
}
