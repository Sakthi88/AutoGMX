#!/usr/bin/env python3
import json
import os
import re
import signal
import shutil
import subprocess
import sys
import uuid
import zipfile
from datetime import datetime
from pathlib import Path

from flask import Flask, abort, flash, redirect, render_template, request, send_file, url_for
from werkzeug.exceptions import RequestEntityTooLarge

# NOTE: Full file restore in progress - if you see this incomplete message, restore from commit 4f292cc and apply the two-line gromos_skeleton patch documented in PR #2.
APP_DIR = Path(__file__).resolve().parent
print('INCOMPLETE_RESTORE')
