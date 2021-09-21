#!/bin/bash

ps -ef | grep covigator-dashboard | grep -v grep | awk '{print $2}' | xargs kill -9
